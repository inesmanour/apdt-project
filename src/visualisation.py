# visualisation.py
"""
Utilitaires de visualisation et d’export PyMOL pour la membrane.

Fonctions fournies
------------------
- _norm : normalisation sûre d’un vecteur.
- _orthobasis_from_normal : construit une base orthonormée (n̂, u, v) à partir
  d’une normale.
- slab_mask : renvoie un masque booléen des points contenus dans la « dalle »
  d’épaisseur 2×half_thick autour du plan normal à n̂ au niveau du COM.
- plot_3d_scene : affiche la scène 3D (Cα, COM, directions candidates,
  meilleure normale et deux dalles) avec Matplotlib.
- _pymol_resi_selection : fabrique une sélection PyMOL robuste à partir d’une
  liste (chaine, resseq, icode).
- write_pml : écrit un script .pml qui charge la protéine, place le COM,
  trace la normale et dessine deux dalles (CGO).
- residues_in_slab_for_pymol : renvoie les résidus tombant dans la dalle
  pour les mettre en évidence dans PyMOL.

Notes
-----
- Aucune dépendance PyMOL pendant l’exécution Python : seul un fichier .pml
  est généré, à charger ensuite dans PyMOL avec `@view_membrane.pml`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence, Tuple, List

import numpy as np


# ---------------------------------------------------------------------
# Outils géométriques
# ---------------------------------------------------------------------
def _norm(v: np.ndarray) -> np.ndarray:
    """Normalise un vecteur.

    Parameters
    ----------
    v : np.ndarray
        Vecteur 3D (ou broadcastable vers (3,)).

    Returns
    -------
    np.ndarray
        Vecteur normalisé. Si la norme vaut 0, retourne le vecteur d’origine.
    """
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v if n == 0 else v / n


def _orthobasis_from_normal(n: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construit une base orthonormée (n̂, u, v) à partir d’une normale.

    L’axe auxiliaire est choisi pour éviter la colinéarité avec la normale.

    Parameters
    ----------
    n : np.ndarray
        Normale (3,) non nécessairement normalisée.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        (n̂, u, v) : n̂ = n/‖n‖, u ⟂ n̂, v = n̂ × u ; tous de norme 1.
    """
    n = _norm(n)
    tmp = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = _norm(np.cross(n, tmp))
    v = _norm(np.cross(n, u))
    return n, u, v


def slab_mask(
    points: np.ndarray,
    com: np.ndarray,
    normal: np.ndarray,
    half_thick: float,
) -> np.ndarray:
    """Sélectionne les points à l’intérieur de la « dalle » centrée sur le plan.

    Un point P est gardé si |(P - COM)·n̂| ≤ half_thick.

    Parameters
    ----------
    points : np.ndarray
        Tableau (N,3) de points 3D.
    com : np.ndarray
        Centre (COM), shape (3,).
    normal : np.ndarray
        Normale au plan (3,). Sera normalisée en interne.
    half_thick : float
        Demi-épaisseur de la dalle (Å).

    Returns
    -------
    np.ndarray
        Masque booléen de taille (N,) : True si le point est dans la dalle.
    """
    n = _norm(normal)
    d = (np.asarray(points, float) - np.asarray(com, float)) @ n
    return np.abs(d) <= float(half_thick)


# ---------------------------------------------------------------------
# Matplotlib (utile pour debug rapide)
# ---------------------------------------------------------------------
def plot_3d_scene(
    coords: np.ndarray,
    com: np.ndarray,
    dirs: Optional[np.ndarray] = None,
    best_normal: Optional[np.ndarray] = None,
    half_thickness: float = 15.0,
    plane_size: float = 60.0,
    highlight_slab: bool = True,
    title: str = "Vue 3D (Cα, COM, directions)",
) -> None:
    """Affiche une scène 3D : Cα, COM, directions (optionnel) et dalles (optionnel).

    Parameters
    ----------
    coords : np.ndarray
        Coordonnées Cα (N,3).
    com : np.ndarray
        Centre de masse (3,).
    dirs : np.ndarray, optional
        Directions candidates (K,3). Tracées depuis COM si fourni.
    best_normal : np.ndarray, optional
        Normale retenue. Si fournie, la normale épaisse + deux dalles sont tracées.
    half_thickness : float, default=15.0
        Demi-épaisseur des dalles (Å).
    plane_size : float, default=60.0
        Taille de l’arête des quadrilatères représentant les dalles (Å).
    highlight_slab : bool, default=True
        Met en évidence les Cα contenus dans la dalle (points dorés).
    title : str, default="Vue 3D (Cα, COM, directions)"
        Titre de la figure.

    Returns
    -------
    None
        Affiche une fenêtre matplotlib.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    coords = np.asarray(coords, float)
    com = np.asarray(com, float)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=8, alpha=0.55, label="Cα")
    ax.scatter([com[0]], [com[1]], [com[2]], s=120, c="r", label="COM")

    extra = [coords, com.reshape(1, 3)]
    if dirs is not None:
        L = np.linalg.norm(coords - com, axis=1).max() * 0.45
        for n in np.asarray(dirs, float):
            n = _norm(n)
            P = com + n * L
            ax.plot([com[0], P[0]], [com[1], P[1]], [com[2], P[2]], lw=1, alpha=0.8, color="#555")
            extra.append(P.reshape(1, 3))

    if best_normal is not None:
        n, u, v = _orthobasis_from_normal(best_normal)
        L = np.linalg.norm(coords - com, axis=1).max() * 0.55
        tip = com + n * L
        ax.plot([com[0], tip[0]], [com[1], tip[1]], [com[2], tip[2]], lw=3, color="#1f77b4", label="normale*")
        extra.append(tip.reshape(1, 3))

        s = plane_size / 2.0
        base = np.stack([-u * s - v * s, u * s - v * s, u * s + v * s, -u * s + v * s], axis=0)
        for sign, fc, ec in [
            (-1, (0.2, 0.6, 1.0, 0.25), (0.2, 0.6, 1.0, 0.8)),
            (1, (1.0, 0.3, 0.3, 0.25), (1.0, 0.3, 0.3, 0.8)),
        ]:
            center = com + sign * n * half_thickness
            quad = center + base
            ax.add_collection3d(Poly3DCollection([quad], facecolors=[fc], edgecolors=[ec], linewidths=1.5))
            extra.append(quad)

        if highlight_slab:
            mask = slab_mask(coords, com, n, half_thickness)
            sel = coords[mask]
            if sel.size:
                ax.scatter(
                    sel[:, 0], sel[:, 1], sel[:, 2],
                    s=32, edgecolor="k", facecolor="gold",
                    alpha=0.95, label="dans la dalle",
                )
                extra.append(sel)

    scene = np.vstack([np.asarray(e).reshape(-1, 3) for e in extra])
    mn, mx = scene.min(0), scene.max(0)
    ax.set_xlim(mn[0], mx[0])
    ax.set_ylim(mn[1], mx[1])
    ax.set_zlim(mn[2], mx[2])
    ax.set_title(title)
    ax.legend(loc="best")
    ax.set_box_aspect((1, 1, 1))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------
# PyMOL
# ---------------------------------------------------------------------
def _pymol_resi_selection(
    triplets: Sequence[Tuple[str, int, str]],
    name_filter: str = "CA",
) -> str:
    """Construit une chaîne de sélection PyMOL robuste.

    Parameters
    ----------
    triplets : Sequence[Tuple[str, int, str]]
        Séquence de tuples (chain, resseq, icode). Utiliser " " pour icode vide.
    name_filter : str, default="CA"
        Nom d’atome à filtrer (ex. "CA").

    Returns
    -------
    str
        Sélection PyMOL (ex. "(chain A and (resi 12+13)) and name CA").
        Retourne "none" si la liste est vide.
    """
    by_chain: dict[str, List[str]] = {}
    for ch, resseq, icode in triplets:
        tok = f"{resseq}" if not icode or icode == " " else f"{resseq}^{icode}"
        by_chain.setdefault(ch, []).append(tok)
    parts = [f"(chain {ch} and (resi {'+'.join(t)}))" for ch, t in by_chain.items()]
    if not parts:
        return "none"
    return f"({' or '.join(parts)}) and name {name_filter}"


def write_pml(
    pdb_path: str | Path,
    com: np.ndarray,
    normal: np.ndarray,
    half_thickness: float = 15.0,
    plane_size: float = 60.0,
    out_pml: str | Path = "view_membrane.pml",
    obj_name: str = "prot",
    slab_residues: Optional[Sequence[Tuple[str, int, str]]] = None,
) -> None:
    """Écrit un script PyMOL (.pml) pour visualiser la membrane estimée.

    Le script :
      - charge le PDB et l’affiche en cartoon,
      - ajoute un pseudo-atome au COM,
      - trace la normale (ligne + cône),
      - dessine deux dalles (polygones CGO) à ± half_thickness le long de la normale,
      - optionnellement sélectionne les Cα dans la dalle et les colore en jaune.

    Parameters
    ----------
    pdb_path : str | Path
        Chemin vers le fichier .pdb à charger dans PyMOL.
    com : np.ndarray
        Centre de masse (3,).
    normal : np.ndarray
        Normale estimée (3,) ; sera normalisée.
    half_thickness : float, default=15.0
        Demi-épaisseur de la dalle (Å).
    plane_size : float, default=60.0
        Taille de l’arête des dalles (Å).
    out_pml : str | Path, default="view_membrane.pml"
        Nom du fichier .pml à écrire.
    obj_name : str, default="prot"
        Nom d’objet PyMOL pour le PDB.
    slab_residues : Sequence[Tuple[str, int, str]] | None, optional
        Résidus (chaine, resseq, icode) à surligner s’ils sont dans la dalle.

    Returns
    -------
    None
        Écrit le fichier `.pml` sur le disque.
    """
    pdb_path = Path(pdb_path)
    out_pml = Path(out_pml)
    com = np.asarray(com, float)
    n, u, v = _orthobasis_from_normal(np.asarray(normal, float))
    s = float(plane_size) / 2.0

    base = np.stack([-u * s - v * s, u * s - v * s, u * s + v * s, -u * s + v * s], axis=0)
    c_in = com - n * half_thickness
    c_out = com + n * half_thickness
    quad_in = c_in + base
    quad_out = c_out + base
    tip = com + n * 25.0

    # triangles (2 par dalle) pour CGO
    tri_in = np.array(
        [quad_in[0], quad_in[1], quad_in[2], quad_in[0], quad_in[2], quad_in[3]],
        dtype=float,
    ).tolist()
    tri_out = np.array(
        [quad_out[0], quad_out[1], quad_out[2], quad_out[0], quad_out[2], quad_out[3]],
        dtype=float,
    ).tolist()
    quad_in_list = quad_in.tolist()
    quad_out_list = quad_out.tolist()

    # sélection optionnelle (pré-calculée côté Python)
    sel_cmd = ""
    if slab_residues:
        sel = _pymol_resi_selection(slab_residues, name_filter="CA")
        sel_cmd = (
            f"select slab_sel, {sel}\n"
            "show spheres, slab_sel\n"
            "set sphere_scale, 0.35, slab_sel\n"
            "color yellow, slab_sel\n"
            "group vis, slab_sel\n"
        )

    def _fmt_vec_list(lst: List[List[float]]) -> str:
        """Formate [[x,y,z], ...] en texte avec 3 décimales (fichier plus léger)."""
        chunks = []
        for x, y, z in lst:
            chunks.append(f"[{x:.3f},{y:.3f},{z:.3f}]")
        return "[" + ",".join(chunks) + "]"

    with out_pml.open("w", encoding="utf-8") as f:
        # commandes de base
        f.write("reinitialize\n")
        f.write(f"load {pdb_path.as_posix()}, {obj_name}\n")
        f.write(f"as cartoon, {obj_name}\n")
        f.write("bg_color white\n")
        f.write("set two_sided_lighting, on\n")
        f.write("set antialias, 2\n")
        f.write("set cgo_lighting, 1\n")
        f.write("set transparency_mode, 2\n")
        f.write("group vis\n")

        # COM
        f.write("pseudoatom com, pos=[%.3f,%.3f,%.3f]\n" % (com[0], com[1], com[2]))
        f.write("show spheres, com\n")
        f.write("set sphere_scale, 0.8, com\n")
        f.write("color red, com\n")
        f.write("group vis, com\n")

        # coins des dalles (repères)
        for i, (x, y, z) in enumerate(quad_in_list):
            f.write(f"pseudoatom pin_{i}, pos=[{x:.3f},{y:.3f},{z:.3f}]\ncolor blue, pin_{i}\n")
            f.write("group vis, pin_%d\n" % i)
        for i, (x, y, z) in enumerate(quad_out_list):
            f.write(f"pseudoatom pout_{i}, pos=[{x:.3f},{y:.3f},{z:.3f}]\ncolor red, pout_{i}\n")
            f.write("group vis, pout_%d\n" % i)

        # CGO (bloc Python)
        f.write("python\n")
        f.write("from pymol import cmd\nfrom pymol.cgo import *\n")
        f.write(f"com = [{com[0]:.3f},{com[1]:.3f},{com[2]:.3f}]\n")
        f.write(f"tip = [{tip[0]:.3f},{tip[1]:.3f},{tip[2]:.3f}]\n")
        f.write(f"tri_in  = {_fmt_vec_list(tri_in)}\n")
        f.write(f"tri_out = {_fmt_vec_list(tri_out)}\n")
        f.write(f"quad_in  = {_fmt_vec_list(quad_in_list)}\n")
        f.write(f"quad_out = {_fmt_vec_list(quad_out_list)}\n")

        f.write(
            "normal_vec = [BEGIN, LINES, COLOR, 0.10,0.45,0.85,"
            f" VERTEX, {com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f},"
            f" VERTEX, {tip[0]:.3f}, {tip[1]:.3f}, {tip[2]:.3f}, END,\n"
            f" CONE, {tip[0]:.3f}, {tip[1]:.3f}, {tip[2]:.3f},"
            f" {tip[0]+n[0]*5:.3f}, {tip[1]+n[1]*5:.3f}, {tip[2]+n[2]*5:.3f},"
            " 1.8, 0.0, 0.10,0.45,0.85, 0.10,0.45,0.85, 1.0, 1.0]\n"
            "cmd.load_cgo(normal_vec, 'normal_vec')\n"
        )

        f.write(
            "def _emit_tri(name, T, rgb):\n"
            "    r,g,b = rgb\n"
            "    obj=[BEGIN, TRIANGLES, COLOR, r,g,b, ALPHA, 0.28]\n"
            "    for x,y,z in T: obj += [VERTEX, float(x),float(y),float(z)]\n"
            "    obj += [END]\n"
            "    cmd.load_cgo(obj, name)\n"
            "def _emit_outline(name, Q, rgb):\n"
            "    r,g,b = rgb\n"
            "    order=[0,1,2,3,0]\n"
            "    obj=[BEGIN, LINES, COLOR, r,g,b]\n"
            "    for i in range(4):\n"
            "        a,b = order[i], order[i+1]\n"
            "        x1,y1,z1 = Q[a]; x2,y2,z2 = Q[b]\n"
            "        obj += [VERTEX, float(x1),float(y1),float(z1), VERTEX, float(x2),float(y2),float(z2)]\n"
            "    obj += [END]\n"
            "    cmd.load_cgo(obj, name)\n"
            "_emit_tri('plane_in',  tri_in,  (0.20,0.60,1.00))\n"
            "_emit_tri('plane_out', tri_out, (1.00,0.30,0.30))\n"
            "_emit_outline('plane_in_outline',  quad_in,  (0.20,0.60,1.00))\n"
            "_emit_outline('plane_out_outline', quad_out, (1.00,0.30,0.30))\n"
        )
        f.write("python end\n")

        # sélection optionnelle
        if sel_cmd:
            f.write(sel_cmd)

        # activation & zoom
        f.write("enable vis\n")
        f.write("enable plane_in\nenable plane_out\nenable plane_in_outline\nenable plane_out_outline\nenable normal_vec\n")
        f.write("zoom vis, buffer=10\n")

    print(f"[PyMOL] Script écrit → {out_pml}  | Dans PyMOL : @{out_pml.name}")


def residues_in_slab_for_pymol(
    residues: Sequence,
    com: np.ndarray,
    normal: np.ndarray,
    half_thickness: float,
) -> List[Tuple[str, int, str]]:
    """Renvoie les résidus dont le Cα tombe dans la dalle.

    Parameters
    ----------
    residues : Sequence
        Objets résidus possédant :
        - coord : np.ndarray (3,)
        - chain : str
        - resid : tuple (hetflag, resseq, icode)
    com : np.ndarray
        Centre de masse (3,).
    normal : np.ndarray
        Normale (3,) ; longueur quelconque.
    half_thickness : float
        Demi-épaisseur de la dalle (Å).

    Returns
    -------
    list[tuple[str, int, str]]
        Triplets (chain, resseq, icode) des résidus dans la dalle.
        L’icode vide est renvoyé comme " ".
    """
    pts = np.stack([r.coord for r in residues], axis=0)
    mask = slab_mask(pts, com, normal, half_thickness)
    keep: List[Tuple[str, int, str]] = []
    for r, m in zip(residues, mask):
        if m:
            het, resseq, icode = r.resid
            keep.append((r.chain, int(resseq), icode if icode else " "))
    return keep