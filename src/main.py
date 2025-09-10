# main.py
"""
Détection (orientation + épaisseur) de membrane pour une protéine donnée.

Pipeline :
1) Lecture des Cα depuis un PDB.
2) Annotation DSSP (ASA absolue) puis calcul RSA normalisée.
3) Échantillonnage de directions (fibonacci sphere).
4) Pour chaque direction : profil le long de la normale, optimisation
   du centre et de l'épaisseur d’une dalle « membrane » en maximisant
   un score hydrophobe pondéré par (1 − RSA), avec garde-fous.
5) Récapitulatif (Q par direction, Q_max, normale, d_in/d_out, centres 3D).
6) Visualisation Matplotlib et export PyMOL (.pml).
7) Enregistre un CSV dans results/ qui reflète exactement les infos imprimées.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import warnings
import csv
import numpy as np

# --- Matplotlib en mode headless pour éviter les PNG blancs ---
import matplotlib
matplotlib.use("Agg")  # backend hors-écran
import matplotlib.pyplot as plt

# --- Réduction d’un warning Biopython (non critique pour ASA/RSA) ---
warnings.filterwarnings("ignore", category=UserWarning, module="Bio.PDB.DSSP")

# --- Imports internes (projet) ---
from protein import charger_calpha, annoter_dssp_text, MAP3TO1
from geometrie import (
    centre_de_masse,
    fibonacci_sphere,
    meilleur_segment_sur_direction,
)
from visualisation import plot_3d_scene, write_pml

# --------------------------------------------------------------------
#  Kyte–Doolittle
# --------------------------------------------------------------------
_KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4,
    "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5, "B": -3.5, "Z": -3.5, "X": 0.0
}

# --------------------------------------------------------------------
#  Helpers E/S
# --------------------------------------------------------------------
def _ensure_dir(p: Path) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p

def _write_console_csv(csv_path: Path, header_rows: list[tuple[str, str]], dir_rows: list[dict]) -> None:
    """Écrit un CSV avec une section SUMMARY (clé/valeur) puis une table DIRECTIONS."""
    csv_path = Path(csv_path)
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        # Résumé
        w.writerow(["SECTION", "KEY", "VALUE"])
        for k, v in header_rows:
            w.writerow(["SUMMARY", k, v])
        w.writerow([])

        # Table directions
        cols = ["rank", "Q", "normal_x", "normal_y", "normal_z", "center", "half", "d_in", "d_out", "is_best"]
        w.writerow(["SECTION"] + cols)
        for row in dir_rows:
            w.writerow(["DIRECTIONS"] + [row[c] for c in cols])

def _save_current_or_returned_fig(fig_or_ax_or_none, out_path: Path) -> None:
    """
    Sauvegarde robuste :
    - si plot_3d_scene retourne (fig, ax) → on prend fig
    - si retourne fig → on prend fig
    - sinon → on prend plt.gcf()
    Puis savefig + close sur CETTE figure (évite PNG blancs).
    """
    fig = None
    obj = fig_or_ax_or_none
    if isinstance(obj, tuple) and len(obj) >= 1 and hasattr(obj[0], "savefig"):
        fig = obj[0]
    elif hasattr(obj, "savefig"):
        fig = obj
    else:
        fig = plt.gcf()

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

# --------------------------------------------------------------------
#  Score hydrophobe « soft » : KD pondéré par (1 − RSA) + garde-fous
# --------------------------------------------------------------------
def make_score_fn(
    thresh: float = 0.20,          # filtre RSA/(1-RSA) mini
    use_inverted_rsa: bool = True, # pondérer par (1 - RSA) (favorise l’intérieur membranaire)
    min_res: int = 12,             # au moins N résidus dans la fenêtre
    min_span: float = 18.0,        # étendue spatiale mini (Å) ~ évite les fenêtres trop fines
    span_penalty: float = 1.0,     # intensité de la pénalisation si span < min_span
):
    """
    score = mean( KD(aa) * weight )  avec weight = (1 − RSA) si use_inverted_rsa=True (sinon RSA)
    Garde-fous :
      - rejette les fenêtres avec < min_res résidus
      - pénalise les fenêtres trop « fines » spatialement (span < min_span)
    """
    def _score(subres):
        if len(subres) < min_res:
            return -1e9

        num = 0.0
        den = 0.0
        coords = []
        for r in subres:
            aa1 = MAP3TO1.get((r.resname3 or "").upper(), "X")
            kd = _KD.get(aa1, 0.0)
            rr = 0.0 if r.rsa is None else float(r.rsa)
            rr = min(1.0, max(0.0, rr))  # clamp

            w = (1.0 - rr) if use_inverted_rsa else rr
            if w >= thresh:
                num += kd * w
                den += w

            if getattr(r, "coord", None) is not None:
                coords.append(r.coord)

        if den <= 1e-12:
            return -1e9

        base = num / den

        # pénalisation si fenêtre trop « courte » spatialement
        if coords:
            arr = np.asarray(coords, dtype=float)
            span = float(np.max(np.ptp(arr, axis=0)))  # max des étendues X/Y/Z
            if span < min_span:
                penalty = span_penalty * (1.0 - (span / min_span))
                base -= penalty

        return base

    return _score


def traiter_un_pdb(
    pdb_path: Path,
    chain: str | None,
    mkdssp: str | None,
    dssp_path: Path | None,
    n_dirs: int,
    pas: float,
    demi_epaisseur: float,
    do_plot: bool,
    pml_out_flag: str | None,
    head: int,
    show_dirs_first: bool,
    widen_step: float,
    max_half_extra: float,
    rsa_thresh: float,
    min_half: float | None,
    max_half: float | None,
    results_dir: Path,
) -> None:
    """Exécute le pipeline complet, imprime le récapitulatif et écrit CSV/PNG/PML dans results/."""
    results_dir = _ensure_dir(results_dir)
    pdbid = Path(pdb_path).stem.upper()
    chain_suffix = f"_{chain}" if chain else ""
    base = f"{pdbid}{chain_suffix}"
    csv_path = results_dir / f"{base}.csv"
    png_path = results_dir / f"{base}_plot.png"
    pml_path = results_dir / f"{base}.pml"
    pre_png_path = results_dir / f"{base}_dirs.png"

    print(f"\n=== {pdb_path.name} ===")

    # 1) Lecture des résidus Cα
    residues = charger_calpha(pdb_path, chain_filter=chain)
    print(f"- Cα lus : {len(residues)}")

    # 2) Annotation DSSP → RSA
    try:
        annoter_dssp_text(pdb_path, residues, dssp_path=dssp_path, mkdssp=mkdssp)
        n_asa = sum(r.asa is not None for r in residues)
        print(f"- ASA (DSSP) : {n_asa}/{len(residues)}")
        if n_asa:
            rs = [r.rsa for r in residues if r.rsa is not None]
            if rs:
                print(f"- RSA moyenne : {np.mean(rs):.3f}")
                if head > 0:
                    print("\n--- Premiers résidus (ASA/RSA) ---")
                    for r in residues[:head]:
                        asa = "NA" if r.asa is None else f"{r.asa:.1f}"
                        rsa = "NA" if r.rsa is None else f"{r.rsa:.3f}"
                        print(f"{r.chain} {r.resid} {r.resname3:3s} | ASA={asa} | RSA={rsa}")
    except Exception as e:
        print(f"[WARN] DSSP indisponible ou erreur : {e}")

    # 3) Centre de masse
    coords = np.stack([r.coord for r in residues], axis=0)
    com = centre_de_masse(r.coord for r in residues)
    print(f"- Centre de masse (Cα) : {com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}")

    # 4) Directions
    dirs = fibonacci_sphere(n_dirs)
    if do_plot and show_dirs_first:
        fig_or_ax = plot_3d_scene(
            coords=coords,
            com=com,
            dirs=dirs,
            best_normal=None,
            half_thickness=demi_epaisseur,
            title=f"{pdb_path.name} — directions initiales",
        )
        try:
            _save_current_or_returned_fig(fig_or_ax, pre_png_path)
            print(f"[OK] PNG (dirs): {pre_png_path}")
        except Exception as e:
            print(f"[WARN] Impossible d'enregistrer la vue 'dirs': {e}")

    # 5) Score
    score_fn = make_score_fn(thresh=rsa_thresh, use_inverted_rsa=True)

    def _as_float(x) -> float:
        arr = np.asarray(x, dtype=float).ravel()
        return float(arr[0]) if arr.size else float(x)

    best_global = (-1e9, None, None, None)  # (Q, n_unit, centre, demi_epaisseur)
    dir_scores: list[tuple[float, np.ndarray, float, float, float, float]] = []

    for normal in dirs:
        ret = meilleur_segment_sur_direction(
            residues=residues,
            coords=coords,
            origin=com,
            normal=normal,
            tranche=1.0,  # balayage des centres à 1 Å (ta consigne)
            pas=pas,
            demi_epaisseur_init=demi_epaisseur,
            score_fn=score_fn,
            widen_step=widen_step,
            max_half_extra=max_half_extra,
        )
        if len(ret) < 3:
            raise RuntimeError("meilleur_segment_sur_direction a renvoyé un tuple inattendu.")
        score, centre, w = ret[:3]

        centre = _as_float(centre)
        w = abs(_as_float(w))

        # Bornes réalistes (si données)
        if (min_half is not None) and (max_half is not None) and (max_half > min_half >= 0.0):
            w = float(np.clip(w, min_half, max_half))

        n_u = np.asarray(normal, float)
        n_u = n_u / (np.linalg.norm(n_u) + 1e-12)

        d_in = centre - w
        d_out = centre + w
        dir_scores.append((float(score), n_u.copy(), centre, w, d_in, d_out))

        if score > best_global[0]:
            best_global = (float(score), n_u.copy(), centre, w)

    # 5bis) Récapitulatif des directions triées par Q décroissant
    dir_scores.sort(key=lambda x: x[0], reverse=True)
    print("\n--- Scores par direction (Q) ---")
    for k, (Q, n_u, c, w, d_in, d_out) in enumerate(dir_scores, start=1):
        n_str = np.array2string(n_u, precision=3, suppress_small=True)
        print(f"{k:2d}. Q={Q:.3f} | n={n_str} | centre={c:.2f} Å | demi-ép={w:.2f} Å | d_in={d_in:.2f} Å | d_out={d_out:.2f} Å")

    # 6) Résultat final
    score_max, n_star, c_star, w_star = best_global
    if n_star is None:
        print("[WARN] Aucune direction n'a produit de score.")
        return

    print(f"\n→ Q_max = {score_max:.3f}")

    def _orthobasis_from_normal(n: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = np.asarray(n, float)
        n = n / (np.linalg.norm(n) + 1e-12)
        tmp = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(n, tmp); u /= (np.linalg.norm(u) + 1e-12)
        v = np.cross(n, u);   v /= (np.linalg.norm(v) + 1e-12)
        return n, u, v

    n_hat, u, v = _orthobasis_from_normal(n_star)
    d_in = c_star - w_star
    d_out = c_star + w_star
    center_in = com + d_in * n_hat
    center_out = com + d_out * n_hat

    print(f"- Normale estimée              : {n_hat}")
    print(f"- Centre optimal (le long de n): {c_star:.2f} Å")
    print(f"- Score max (KD×(1−RSA) moy.)  : {score_max:.3f}")
    print(f"- Demi-épaisseur optimale      : {w_star:.2f} Å (épaisseur ≈ {2*w_star:.2f} Å)")
    print(f"- Positions d_in / d_out       : {d_in:.2f} Å / {d_out:.2f} Å le long de n")
    print(f"- Centres des dalles (3D)      : IN  = [{center_in[0]:.3f}, {center_in[1]:.3f}, {center_in[2]:.3f}] | OUT = [{center_out[0]:.3f}, {center_out[1]:.3f}, {center_out[2]:.3f}]")

    # 7) Visualisation finale (+ sauvegarde PNG dans results/)
    if do_plot:
        fig_or_ax = plot_3d_scene(
            coords=coords,
            com=com,
            dirs=None,
            best_normal=n_star,
            half_thickness=w_star,
            title=f"{pdb_path.name} — normale & dalles",
        )
        try:
            _save_current_or_returned_fig(fig_or_ax, png_path)
            print(f"[OK] PNG: {png_path}")
        except Exception as e:
            print(f"[WARN] Impossible d'enregistrer le PNG: {e}")

    # 8) Export PyMOL dans results/<PDB>[_<chain>].pml
    if pml_out_flag:
        write_pml(
            pdb_path=pdb_path,
            com=com,
            normal=n_hat,
            half_thickness=w_star,
            plane_size=60.0,
            out_pml=pml_path,
        )
        print(f"[OK] PML: {pml_path}")

    # 9) CSV miroir de la console dans results/<PDB>[_<chain>].csv
    header_rows = [
        ("pdb", pdbid),
        ("chain", str(chain) if chain else "None"),
        ("n_ca", str(len(residues))),
        ("rsa_thresh", f"{rsa_thresh:.3f}"),
        ("n_dirs", str(n_dirs)),
        ("pas_centres_A", "1.0"),  # balayage des centres à 1 Å
        ("widen_step_A", f"{widen_step:.3f}"),
        ("max_half_extra_A", f"{max_half_extra:.3f}"),
        ("min_half_A", "None" if min_half is None else f"{min_half:.3f}"),
        ("max_half_A", "None" if max_half is None else f"{max_half:.3f}"),
        ("COM_x", f"{com[0]:.3f}"),
        ("COM_y", f"{com[1]:.3f}"),
        ("COM_z", f"{com[2]:.3f}"),
        ("Q_max", f"{score_max:.3f}"),
        ("best_center_A", f"{c_star:.2f}"),
        ("best_half_A", f"{w_star:.2f}"),
        ("best_thickness_A", f"{2*w_star:.2f}"),
        ("best_normal_x", f"{n_hat[0]:.6f}"),
        ("best_normal_y", f"{n_hat[1]:.6f}"),
        ("best_normal_z", f"{n_hat[2]:.6f}"),
        ("d_in_A", f"{d_in:.2f}"),
        ("d_out_A", f"{d_out:.2f}"),
        ("center_in_x", f"{center_in[0]:.3f}"),
        ("center_in_y", f"{center_in[1]:.3f}"),
        ("center_in_z", f"{center_in[2]:.3f}"),
        ("center_out_x", f"{center_out[0]:.3f}"),
        ("center_out_y", f"{center_out[1]:.3f}"),
        ("center_out_z", f"{center_out[2]:.3f}"),
    ]
    dir_rows = []
    for k, (Q, n_u, c, w, din, dout) in enumerate(dir_scores, start=1):
        dir_rows.append({
            "rank": k,
            "Q": round(Q, 6),
            "normal_x": round(float(n_u[0]), 6),
            "normal_y": round(float(n_u[1]), 6),
            "normal_z": round(float(n_u[2]), 6),
            "center": round(float(c), 3),
            "half": round(float(w), 3),
            "d_in": round(float(din), 3),
            "d_out": round(float(dout), 3),
            "is_best": 1 if k == 1 else 0,
        })
    _write_console_csv(csv_path, header_rows, dir_rows)
    print(f"[OK] CSV: {csv_path}")


def main() -> None:
    """Parse les arguments CLI, puis lance le pipeline."""
    ap = argparse.ArgumentParser(
        description="Détection TM (DSSP ACC→RSA, tranches 1 Å, optimisation de fenêtre)"
    )
    ap.add_argument("--pdb", type=Path, help="Un fichier .pdb à traiter")
    ap.add_argument("--dssp", type=Path, default=None,
                    help="Fichier .dssp déjà généré (optionnel). Sinon mkdssp sera appelé.")
    ap.add_argument("--membrane-dir", type=Path, default=Path("data/membrane"),
                    help="Dossier de PDB membranaires (batch si --pdb absent).")
    ap.add_argument("--globular-dir", type=Path, default=Path("data/globular"),
                    help="Dossier de PDB globulaires (batch si --pdb absent).")
    ap.add_argument("--chain", type=str, default=None, help="Chaîne spécifique (ex: A)")
    ap.add_argument("--mkdssp", type=str, default=None,
                    help="Chemin explicite vers mkdssp/dssp (sinon auto-détection).")
    ap.add_argument("--n_dirs", type=int, default=64, help="Nombre de directions sur la sphère")
    ap.add_argument("--pas", type=float, default=1.0, help="Pas de déplacement de la dalle (Å)")
    ap.add_argument("--demi_epaisseur", type=float, default=15.0,
                    help="Demi-épaisseur initiale (Å) de la dalle (point de départ).")
    ap.add_argument("--plot", action="store_true", help="Affiche la scène 3D Matplotlib")
    ap.add_argument("--dirs-first", action="store_true",
                    help="Si --plot, sauvegarde aussi la vue COM+directions avant scoring.")
    ap.add_argument("--pml", type=str, default=None,
                    help="Activer l'export PyMOL (le fichier sera nommé automatiquement).")
    ap.add_argument("--head", type=int, default=0,
                    help="Affiche les N premiers résidus (ASA/RSA) pour contrôle.")
    ap.add_argument("--widen-step", type=float, default=0.5,
                    help="Pas d'élargissement de la membrane (Å) lors de l’optimisation.")
    ap.add_argument("--max-half-extra", type=float, default=2.0,
                    help="Demi-épaisseur additionnelle max (Å) autorisée lors de l'élargissement.")
    # Garde-fous
    ap.add_argument("--rsa-thresh", type=float, default=0.20,
                    help="Seuil minimal (RSA ou 1−RSA) pour contribuer au score.")
    ap.add_argument("--no-widen", action="store_true",
                    help="Désactive l'élargissement (équivaut à --max-half-extra 0).")
    ap.add_argument("--min-half", type=float, default=10.0,
                    help="Borne basse (Å) de la demi-épaisseur (None: -1).")
    ap.add_argument("--max-half", type=float, default=20.0,
                    help="Borne haute (Å) de la demi-épaisseur (None: -1).")
    # Sorties
    ap.add_argument("--results-dir", type=Path, default=Path("results"),
                    help="Dossier de sortie (CSV/PNG/PML)")

    args = ap.parse_args()

    if args.no_widen:
        args.max_half_extra = 0.0

    # Désactiver bornes via -1
    min_half = None if (args.min_half is not None and args.min_half < 0) else args.min_half
    max_half = None if (args.max_half is not None and args.max_half < 0) else args.max_half

    if args.pdb:
        traiter_un_pdb(
            pdb_path=args.pdb,
            chain=args.chain,
            mkdssp=args.mkdssp,
            dssp_path=args.dssp,
            n_dirs=args.n_dirs,
            pas=args.pas,
            demi_epaisseur=args.demi_epaisseur,
            do_plot=args.plot,
            pml_out_flag=args.pml,
            head=args.head,
            show_dirs_first=args.dirs_first,
            widen_step=args.widen_step,
            max_half_extra=args.max_half_extra,
            rsa_thresh=args.rsa_thresh,
            min_half=min_half,
            max_half=max_half,
            results_dir=args.results_dir,
        )
    else:
        for label, folder in (("membrane", args.membrane_dir), ("globular", args.globular_dir)):
            if not folder.exists():
                continue
            print(f"\n### Groupe {label} ({folder})")
            for pdb in sorted(folder.glob("*.pdb")):
                traiter_un_pdb(
                    pdb_path=pdb,
                    chain=args.chain,
                    mkdssp=args.mkdssp,
                    dssp_path=args.dssp,
                    n_dirs=args.n_dirs,
                    pas=args.pas,
                    demi_epaisseur=args.demi_epaisseur,
                    do_plot=args.plot,
                    pml_out_flag=args.pml,
                    head=args.head,
                    show_dirs_first=args.dirs_first,
                    widen_step=args.widen_step,
                    max_half_extra=args.max_half_extra,
                    rsa_thresh=args.rsa_thresh,
                    min_half=min_half,
                    max_half=max_half,
                    results_dir=args.results_dir,
                )


if __name__ == "__main__":
    main()