# hydrophobicite.py
"""
Scores d’hydrophobicité (Kyte–Doolittle) et pondération par RSA.

Ce module fournit :
- un dictionnaire KD (hydrophobicité par acide aminé, code 1 lettre),
- un mapping 3 lettres → 1 lettre (MAP3TO1),
- deux fonctions utilitaires :
    * score_hydro : retourne le score KD d’un acide aminé,
    * weighted_score : retourne un score KD pondéré par l’exposition RSA.

"""

from __future__ import annotations
from typing import Dict, Optional

__all__ = ["KD", "MAP3TO1", "score_hydro", "weighted_score"]

# ---------------------------------------------------------------------
# Tables de correspondance
# ---------------------------------------------------------------------
KD: Dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

MAP3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


# ---------------------------------------------------------------------
# Fonctions
# ---------------------------------------------------------------------
def score_hydro(resname3: str) -> float:
    """Retourne le score d’hydrophobicité Kyte–Doolittle pour un résidu.

    Parameters
    ----------
    resname3 : str
        Code acide aminé sur 3 lettres (ex. "LEU", "GLY"). La casse est ignorée.

    Returns
    -------
    float
        Score KD en unités Kyte–Doolittle. Retourne 0.0 si le résidu n’est pas reconnu.

    Examples
    --------
    >>> score_hydro("LEU")
    3.8
    >>> score_hydro("ASP")
    -3.5
    """
    aa1 = MAP3TO1.get((resname3 or "").upper())
    return KD.get(aa1, 0.0) if aa1 else 0.0


def weighted_score(resname3: str, rsa: Optional[float], thresh: float = 0.25) -> float:
    """Retourne un score KD pondéré par l’exposition (RSA).

    Un résidu est compté uniquement s’il est suffisamment exposé :
    - si RSA est None ou RSA < `thresh`, score = 0.0 (ignoré),
    - sinon, score = score_hydro(resname3).

    Parameters
    ----------
    resname3 : str
        Code acide aminé 3 lettres (ex. "LEU").
    rsa : float | None
        Relative Solvent Accessibility (0..1). None si indisponible.
    thresh : float, default=0.25
        Seuil RSA en dessous duquel on ignore le résidu.

    Returns
    -------
    float
        0.0 si le résidu est enfoui, sinon le score KD.

    Examples
    --------
    >>> weighted_score("LEU", rsa=0.30, thresh=0.25)
    3.8
    >>> weighted_score("LEU", rsa=0.10, thresh=0.25)
    0.0
    """
    if rsa is None or rsa < thresh:
        return 0.0
    return score_hydro(resname3)