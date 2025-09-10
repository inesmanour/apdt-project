# protein.py
"""
Gestion des résidus Cα et annotation ASA/RSA via DSSP.

Ce module fournit :
- un modèle léger de résidu porté par le Cα (`ResidueCA`),
- des utilitaires pour charger les Cα d'un fichier PDB,
- l’intégration avec `mkdssp/dssp` (exécution et parsing texte),
- le calcul de RSA = ACC / ACC_MAX et l’annotation sur la liste de résidus.

Flux typique
------------
1) `charger_calpha(pdb_path)` extrait les Cα (coordonnées 3D, nom AA, etc.).
2) `annoter_dssp_text(pdb_path, residues)` lit (ou génère) un `.dssp`,
   extrait l’ACC absolu et calcule la RSA, puis remplit `ResidueCA.asa/rsa`.

Notes
-----
- Le module ne dépend pas de DSSP pour *lire* le PDB (uniquement pour
  l'annotation ASA/RSA).
- Le parsing `.dssp` utilise le format texte traditionnel.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import subprocess
import shutil
import tempfile

import numpy as np
from Bio.PDB import PDBParser, is_aa

__all__ = [
    "ACC_MAX",
    "MAP3TO1",
    "ResidueCA",
    "charger_calpha",
    "parse_dssp_text",
    "annoter_dssp_text",
]

# ---------------------------------------------------------------------
# Tables ACC maximal (Tien/Miller) et mappings AA
# ---------------------------------------------------------------------
ACC_MAX: Dict[str, float] = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223,
    "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
    "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174,
}

# 3 lettres → 1 lettre (standard)
MAP3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

# Variantes rencontrées dans les PDB/DSSP → 1 lettre
VARIANT_1LETTER: Dict[str, str] = {
    "HSD": "H", "HSE": "H", "HSP": "H",  # histidines protonations
    "MSE": "M",                          # sélénométhionine
    "ASX": "B",                          # Asn/Asp indéterminé
    "GLX": "Z",                          # Gln/Glu indéterminé
}

# ---------------------------------------------------------------------
# Modèle de résidu (Cα)
# ---------------------------------------------------------------------
@dataclass
class ResidueCA:
    """Représentation minimale d’un résidu porté par son atome Cα.

    Attributes
    ----------
    chain : str
        ID de chaîne (ex. "A").
    resid : tuple[str, int, str]
        Identifiant PDB du résidu : (hetflag, resseq, icode).
        Exemple : (' ', 123, ' ') pour un résidu standard sans insercode.
    resname3 : str
        Nom de l’acide aminé en 3 lettres (ex. "LEU").
    coord : np.ndarray
        Coordonnées 3D (shape (3,), en Å) de l’atome Cα.
    asa : float | None
        Surface accessible (ACC absolu, en Å²) si connue.
    rsa : float | None
        Relative Solvent Accessibility (ACC / ACC_MAX), bornée dans [0, 1]
        si `annoter_dssp_text(..., clamp=True)`.
    """
    chain: str
    resid: Tuple[str, int, str]
    resname3: str
    coord: np.ndarray
    asa: Optional[float] = None
    rsa: Optional[float] = None

# ---------------------------------------------------------------------
# Utilitaires fichiers / exécution
# ---------------------------------------------------------------------
def _assert_file_exists(p: Path) -> None:
    """Vérifie l’existence d’un fichier.

    Parameters
    ----------
    p : Path
        Chemin du fichier à vérifier.

    Raises
    ------
    FileNotFoundError
        Si le fichier n’existe pas.
    """
    if not p.exists():
        raise FileNotFoundError(f"Fichier introuvable: {p}")


def _auto_mkdssp() -> Optional[str]:
    """Recherche `mkdssp` ou `dssp` dans le PATH.

    Returns
    -------
    str | None
        Chemin exécutable si trouvé, sinon `None`.
    """
    for exe in ("mkdssp", "dssp"):
        path = shutil.which(exe)
        if path:
            return path
    return None

# ---------------------------------------------------------------------
# Extraction des Cα depuis un PDB
# ---------------------------------------------------------------------
def charger_calpha(
    pdb_path: str | Path,
    model_index: int = 0,
    chain_filter: str | None = None,
) -> List[ResidueCA]:
    """Charge tous les atomes Cα d’un PDB et construit les `ResidueCA`.

    Parameters
    ----------
    pdb_path : str | Path
        Chemin du fichier `.pdb`.
    model_index : int, default=0
        Index du modèle à lire (structure multi-modèles).
    chain_filter : str | None, default=None
        Si défini, ne garde que cette chaîne (ex. "A").

    Returns
    -------
    list[ResidueCA]
        Liste triée par (chaîne, resseq, icode) des résidus portés par leur Cα.

    Raises
    ------
    FileNotFoundError
        Si le PDB n’existe pas.
    IndexError
        Si `model_index` n’est pas présent dans la structure.
    RuntimeError
        Si aucun Cα n’est trouvé (fichier/chaîne incorrects).
    """
    p = Path(pdb_path)
    _assert_file_exists(p)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", str(p))
    try:
        model = list(structure)[model_index]
    except Exception as e:  # pragma: no cover - message explicite
        raise IndexError(
            f"model_index={model_index} indisponible dans {p.name}"
        ) from e

    out: List[ResidueCA] = []
    chain_norm = chain_filter.upper() if chain_filter else None

    for chain in model:
        if chain_norm and chain.id.upper() != chain_norm:
            continue
        for res in chain:
            # standard=False pour garder MSE, etc.; on ne garde que les résidus AA
            if not is_aa(res, standard=False):
                continue
            if "CA" not in res:
                continue
            out.append(
                ResidueCA(
                    chain=chain.id,
                    resid=res.get_id(),  # (' ', resseq, icode)
                    resname3=res.get_resname(),
                    coord=res["CA"].get_coord().astype(float),
                )
            )

    if not out:
        raise RuntimeError(
            f"Aucun Cα trouvé (pdb={p.name}, chain={chain_filter or 'ALL'})"
        )

    out.sort(key=lambda r: (r.chain, r.resid[1], r.resid[2]))
    return out

# ---------------------------------------------------------------------
# Appel à mkdssp et parsing du .dssp texte
# ---------------------------------------------------------------------
def _run_mkdssp(pdb_path: Path, out_dssp: Path, mkdssp: Optional[str]) -> None:
    """Exécute mkdssp/dssp pour produire un fichier `.dssp`.

    Tente plusieurs syntaxes pour couvrir les différentes versions :
    1) `mkdssp -i IN -o OUT`
    2) `mkdssp IN OUT`
    3) `mkdssp` avec redirections (stdin->stdout)

    Parameters
    ----------
    pdb_path : Path
        Fichier PDB d’entrée.
    out_dssp : Path
        Chemin de sortie du `.dssp` texte.
    mkdssp : str | None
        Chemin explicite vers l’exécutable. Si `None`, recherche auto.

    Raises
    ------
    RuntimeError
        Si aucune syntaxe ne fonctionne ou si l’exécutable est introuvable.
    """
    if mkdssp is None:
        mkdssp = _auto_mkdssp()
    if not mkdssp:
        raise RuntimeError(
            "mkdssp/dssp introuvable. Installez `dssp` et/ou ajoutez-le au PATH."
        )

    # 1) mkdssp -i IN -o OUT
    try:
        subprocess.run(
            [mkdssp, "-i", str(pdb_path), "-o", str(out_dssp)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if out_dssp.exists():
            return
    except subprocess.CalledProcessError:
        pass

    # 2) mkdssp IN OUT
    try:
        subprocess.run(
            [mkdssp, str(pdb_path), str(out_dssp)],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if out_dssp.exists():
            return
    except subprocess.CalledProcessError:
        pass

    # 3) mkdssp (stdin->stdout) via pipes
    with open(pdb_path, "rb") as fin, open(out_dssp, "wb") as fout:
        res = subprocess.run([mkdssp], stdin=fin, stdout=fout, stderr=subprocess.PIPE)
    if res.returncode == 0 and out_dssp.exists():
        return

    stderr = res.stderr.decode(errors="ignore") if res.stderr else ""
    raise RuntimeError(
        f"mkdssp a échoué (code {res.returncode}). stderr:\n{stderr}"
    )


def _aa3_to_aa1(aa3_or_aa1: str) -> str:
    """Convertit un acide aminé 3 lettres (ou déjà 1 lettre) en 1 lettre.

    Gère les variantes les plus courantes (`VARIANT_1LETTER`). Retourne
    `'X'` si la lettre n’est pas reconnue.

    Parameters
    ----------
    aa3_or_aa1 : str
        Code AA potentiel (1 ou 3 lettres).

    Returns
    -------
    str
        Code 1 lettre (A, R, N, …, X).
    """
    a = (aa3_or_aa1 or "").strip().upper()
    if len(a) == 1:
        return a if (a in ACC_MAX or a in {"B", "Z", "U", "O", "X", "!"}) else "X"
    if a in VARIANT_1LETTER:
        return VARIANT_1LETTER[a]
    return MAP3TO1.get(a, "X")


def parse_dssp_text(
    dssp_path: str | Path,
) -> Dict[Tuple[str, Tuple[str, int, str]], Tuple[float, str]]:
    """Parse un fichier `.dssp` **texte** et retourne un index ACC/AA1.

    L’index associe chaque entrée (chaîne, resid) → (ACC_absolu, AA_1lettre).

    Le parser repose sur le format historique DSSP :
    - `resseq` : colonnes 5:10
    - `icode`  : colonne 10
    - `chain`  : colonne 11
    - `aa`     : colonne 13 (souvent 1 lettre, `'!'` = break)
    - `ACC`    : colonnes 34:38 (parfois 34:39)

    Parameters
    ----------
    dssp_path : str | Path
        Chemin du fichier `.dssp` texte.

    Returns
    -------
    dict
        Dictionnaire : (chain, (' ', resseq, icode)) → (ACC_absolu, AA_1lettre)

    Raises
    ------
    FileNotFoundError
        Si le fichier n’existe pas.
    """
    dssp_path = Path(dssp_path)
    _assert_file_exists(dssp_path)

    idx: Dict[Tuple[str, Tuple[str, int, str]], Tuple[float, str]] = {}
    in_table = False

    with dssp_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not in_table:
                if line.startswith("  #  RESIDUE"):
                    in_table = True
                continue
            if len(line) < 40:
                continue

            try:
                resseq_str = line[5:10].strip()
                icode = line[10].strip() or " "
                chain = line[11].strip() or " "
                aa_field = line[13].strip()  # 'A'… ou '!' (break)

                if not resseq_str:
                    continue
                resseq = int(resseq_str)

                acc_str = line[34:38].strip() or line[34:39].strip()
                acc = float(acc_str) if acc_str else 0.0

                aa1 = "X" if aa_field == "!" else _aa3_to_aa1(aa_field)
            except Exception:
                # Ligne mal formée → ignorée
                continue

            hetflag = " "
            res_id = (hetflag, int(resseq), icode)
            idx[(chain, res_id)] = (acc, aa1)

    return idx

# ---------------------------------------------------------------------
# Annotation ASA/RSA sur la liste de Cα
# ---------------------------------------------------------------------
def annoter_dssp_text(
    pdb_path: str | Path,
    residues: List[ResidueCA],
    dssp_path: str | Path | None = None,
    mkdssp: str | None = None,
    clamp: bool = True,
) -> None:
    """Annote `residues` avec ASA (ACC) et RSA à partir d’un `.dssp`.

    Étapes :
    1) Lire le `.dssp` (ou l’engendrer via `mkdssp`) pour obtenir
       l’ACC absolu (Å²) et l’acide aminé en 1 lettre.
    2) Calculer `RSA = ACC / ACC_MAX[AA1]`.
    3) Si `clamp=True`, borner RSA dans l’intervalle [0, 1].
    4) Renseigner `r.asa` et `r.rsa` pour chaque `ResidueCA`.

    Parameters
    ----------
    pdb_path : str | Path
        PDB utilisé par DSSP (sert aussi si `dssp_path` est None).
    residues : list[ResidueCA]
        Résidus Cα à annoter (modifiés en place).
    dssp_path : str | Path | None, default=None
        Chemin d’un fichier `.dssp` **déjà** généré. Si `None`,
        `mkdssp` sera exécuté dans un dossier temporaire.
    mkdssp : str | None, default=None
        Chemin explicite vers l’exécutable. Si `None`, recherche automatique.
    clamp : bool, default=True
        Si vrai, clip RSA dans [0, 1].

    Raises
    ------
    FileNotFoundError
        Si `pdb_path` n’existe pas.
    RuntimeError
        Si `mkdssp` est introuvable/échoue lorsqu’il est nécessaire.
    """
    pdb_path = Path(pdb_path)
    _assert_file_exists(pdb_path)

    # 1) obtenir l’index ACC/AA1
    if dssp_path is None:
        # Générer dans un fichier temporaire puis parser
        with tempfile.TemporaryDirectory() as tmpd:
            out_dssp = Path(tmpd) / (pdb_path.stem + ".dssp")
            _run_mkdssp(pdb_path, out_dssp, mkdssp)
            idx = parse_dssp_text(out_dssp)
    else:
        idx = parse_dssp_text(dssp_path)

    # 2) appliquer aux résidus
    missing = 0
    for r in residues:
        hetflag, resseq, icode = r.resid
        icode = icode if (icode and icode != "") else " "
        key = (r.chain, (hetflag, int(resseq), icode))

        acc_abs, aa1 = idx.get(key, (None, None))
        if acc_abs is None or aa1 is None:
            r.asa = None
            r.rsa = None
            missing += 1
            continue

        max_acc = ACC_MAX.get(aa1)
        rsa = (acc_abs / float(max_acc)) if (max_acc and max_acc > 0) else 0.0
        if clamp:
            rsa = max(0.0, min(1.0, rsa))

        r.asa = float(acc_abs)
        r.rsa = float(rsa)

    if missing:
        print(
            f"[WARN] DSSP: {missing} résidu(s) non apparié(s) "
            "(chaîne/resseq/icode)."
        )