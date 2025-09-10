# geometrie.py
"""
Outils géométriques pour la détection de membrane.

Ce module regroupe :
- des utilitaires de géométrie (centre de masse, normalisation, base
  orthonormée implicite),
- un échantillonnage quasi-uniforme de directions sur la sphère
  (fibonacci_sphere),
- la recherche, **le long d’une direction donnée**, de la dalle (slab)
  la plus “hydrophobe” selon un score fourni par l’utilisateur
  (meilleur_segment_sur_direction).

Le cœur de l’algorithme est *meilleur_segment_sur_direction* :
1) projeter tous les Cα sur l’axe (origin, normal),
2) balayer des centres de fenêtre tous les `tranche` Å pour une
   demi-épaisseur initiale `demi_epaisseur_init`,
3) à partir du meilleur centre, **élargir** la fenêtre par pas
   `widen_step` tant que le score augmente, avec une borne maximale
   `max_half_extra` d’augmentation de demi-épaisseur.

Notes
-----
- Le *score* est fourni par l’appelant via `score_fn(subres)` où
  *subres* est la **sous-liste de résidus** à l’intérieur de la fenêtre.
- La forme des résidus n’est **pas imposée** : on suppose seulement
  qu’ils possèdent au minimum l’attribut `.coord` (np.ndarray(3,)).
- Le score est libre (hydrophobicité, RSA pondéré, etc.).
"""

from __future__ import annotations

from typing import Callable, Iterable, List, Tuple
import numpy as np

__all__ = [
    "centre_de_masse",
    "fibonacci_sphere",
    "meilleur_segment_sur_direction",
]

# ---------------------------------------------------------------------
# Utilitaires
# ---------------------------------------------------------------------
def centre_de_masse(coords_iter: Iterable[np.ndarray]) -> np.ndarray:
    """Calcule le centre de masse (moyenne) d’un ensemble de coordonnées 3D.

    Parameters
    ----------
    coords_iter : Iterable[np.ndarray]
        Itérable de vecteurs 3D (shape (3,)).

    Returns
    -------
    np.ndarray
        Vecteur 3D (shape (3,)) correspondant à la moyenne.
    """
    coords = np.stack(list(coords_iter), axis=0)
    return coords.mean(axis=0)


def fibonacci_sphere(n: int) -> np.ndarray:
    """Génère `n` directions quasi-uniformes sur la sphère unité.

    Implémentation classique via l’angle d’or (Fibonacci lattices).

    Parameters
    ----------
    n : int
        Nombre de directions souhaité (≥ 1).

    Returns
    -------
    np.ndarray
        Tableau (n, 3) de vecteurs **unitaires**.
    """
    pts = []
    phi = (1.0 + 5 ** 0.5) / 2.0
    ga = 2.0 * np.pi * (1.0 - 1.0 / phi)
    for i in range(n):
        z = 1.0 - (2 * i + 1) / n
        r = float(np.sqrt(max(0.0, 1.0 - z * z)))
        th = i * ga
        pts.append([r * np.cos(th), r * np.sin(th), z])
    pts = np.asarray(pts, float)
    pts /= (np.linalg.norm(pts, axis=1, keepdims=True) + 1e-12)
    return pts


def _unit(v: np.ndarray) -> np.ndarray:
    """Retourne le vecteur unitaire colinéaire à `v` (robuste aux très petites normes)."""
    v = np.asarray(v, float)
    return v / (np.linalg.norm(v) + 1e-12)


def _profile_centers(t: np.ndarray, tranche: float) -> np.ndarray:
    """Construit les centres de fenêtres le long de l’axe projeté.

    On étend la plage min–max des projections avec une marge pour
    s’assurer que la dalle peut couvrir toute la protéine.

    Parameters
    ----------
    t : np.ndarray
        Projections scalaires des Cα sur l’axe (origin, normal).
    tranche : float
        Pas entre deux centres successifs (en Å).

    Returns
    -------
    np.ndarray
        Tableau 1D des positions de centres (peut être vide).
    """
    # marge pour couvrir tout l’objet (min 4 Å ou 15 % de l’étendue)
    margin = max(4.0, 0.15 * (t.max() - t.min()))
    tmin = float(np.floor(t.min() - margin))
    tmax = float(np.ceil(t.max() + margin))
    edges = np.arange(tmin, tmax + tranche, tranche, dtype=float)
    if edges.size < 2:
        return np.array([], dtype=float)
    return (edges[:-1] + edges[1:]) / 2.0


def _score_window(
    residues: List,
    t: np.ndarray,
    center: float,
    half: float,
    score_fn: Callable[[List], float],
) -> float:
    """Calcule le score de la *fenêtre* centrée en `center` de demi-épaisseur `half`.

    Parameters
    ----------
    residues : list
        Liste de résidus (doivent au moins posséder `.coord`).
    t : np.ndarray
        Projections des résidus sur l’axe (origin, normal).
    center : float
        Position du centre de la fenêtre (en Å) le long de l’axe.
    half : float
        Demi-épaisseur de la fenêtre (en Å).
    score_fn : Callable[[list], float]
        Fonction utilisateur `score_fn(subres)`.

    Returns
    -------
    float
        Score numérique (−inf si aucun point dans la fenêtre).
    """
    mask = np.abs(t - center) <= half
    if not np.any(mask):
        return -np.inf
    sub = [residues[i] for i, m in enumerate(mask) if m]
    return float(score_fn(sub))


# ---------------------------------------------------------------------
# Cœur : recherche de la dalle optimale sur une direction donnée
# ---------------------------------------------------------------------
def meilleur_segment_sur_direction(
    residues: List,
    coords: np.ndarray,
    origin: np.ndarray,
    normal: np.ndarray,
    tranche: float = 1.0,
    pas: float = 1.0,  # non utilisé ici (compat API), le balayage se fait avec `tranche`
    demi_epaisseur_init: float = 15.0,
    score_fn: Callable[[List], float] = lambda rs: 0.0,
    widen_step: float = 1.0,
    max_half_extra: float = 30.0,
) -> Tuple[float, np.ndarray, float, float]:
    """Trouve la dalle de score maximal **le long d’une direction**.

    Étapes
    ------
    1) Projeter toutes les coordonnées `coords` sur l’axe défini par
       `origin` et `normal`.
    2) Balayer des centres de fenêtre espacés de `tranche` Å, en gardant
       la demi-épaisseur **initiale** `demi_epaisseur_init`. On retient
       le centre donnant le meilleur score.
    3) À partir de ce centre, **élargir** la fenêtre par pas `widen_step`
       tant que le score augmente. L’augmentation totale de demi-épaisseur
       est bornée par `max_half_extra`.

    Paramètres
    ----------
    residues : list
        Liste de résidus (objets avec `.coord` requis ; autres champs libres).
    coords : np.ndarray
        Tableau (N, 3) des coordonnées 3D (souvent Cα). Doit correspondre à `residues`.
    origin : np.ndarray
        Point d’origine (3D) pour définir l’axe de projection.
    normal : np.ndarray
        Vecteur direction (3D) pour l’axe (sera normalisé).
    tranche : float, default=1.0
        Pas entre les **centres** testés (en Å). Contrôle la granularité.
    pas : float, default=1.0
        Paramètre conservé pour compatibilité (le balayage utilise `tranche`).
    demi_epaisseur_init : float, default=15.0
        Demi-épaisseur initiale (en Å) de la fenêtre lors du balayage.
    score_fn : Callable[[list], float], default=lambda rs: 0.0
        Fonction de score appliquée à la sous-liste de résidus à l’intérieur
        de la fenêtre courante.
    widen_step : float, default=1.0
        Incrément d’élargissement de la demi-épaisseur (en Å) pendant l’étape 3.
    max_half_extra : float, default=30.0
        Demi-épaisseur **additionnelle maximale** (en Å) pouvant être ajoutée
        lors de l’élargissement.

    Returns
    -------
    tuple
        `(score_final, n_unitaire, centre_optimal, demi_epaisseur_opt)` où :
        - `score_final` : float, score obtenu sur la fenêtre optimale,
        - `n_unitaire`  : np.ndarray(3,), la direction normalisée,
        - `centre_optimal` : float, position (en Å) le long de `n_unitaire`,
        - `demi_epaisseur_opt` : float, demi-épaisseur finale (en Å).

    Notes
    -----
    - En cas d’absence de centres (projections dégénérées), la fonction
      retourne `(-inf, n, 0.0, demi_epaisseur_init)` pour signaler l’échec.
    """
    # normalise la direction et projette les points
    n = _unit(normal)
    origin = np.asarray(origin, float)
    t = (coords - origin) @ n

    # centres candidats le long de la normale
    centers = _profile_centers(t, tranche=tranche)
    if centers.size == 0:
        return (-np.inf, n, 0.0, float(demi_epaisseur_init))

    # 1) meilleur centre pour la demi-épaisseur initiale
    half = float(demi_epaisseur_init)
    best_center = None
    best_score = -np.inf
    for c in centers:
        s = _score_window(residues, t, c, half, score_fn)
        if s > best_score:
            best_score, best_center = s, float(c)

    if best_center is None:
        return (-np.inf, n, 0.0, half)

    # 2) élargissement autour du meilleur centre
    def score_bounds(d_in: float, d_out: float) -> float:
        c = 0.5 * (d_in + d_out)
        h = 0.5 * (d_out - d_in)
        return _score_window(residues, t, c, h, score_fn)

    c_star = best_center
    d_in = c_star - half
    d_out = c_star + half
    base = best_score

    limit = half + float(max_half_extra)
    improved = True
    eps = 1e-12

    while improved:
        improved = False

        # (a) élargissement symétrique
        if (0.5 * (d_out - d_in)) + widen_step <= limit:
            d_in_try = d_in - widen_step
            d_out_try = d_out + widen_step
            s_try = score_bounds(d_in_try, d_out_try)
            if s_try > base + eps:
                d_in, d_out, base = d_in_try, d_out_try, s_try
                improved = True
                continue

        # (b) déplacer IN seul
        d_in_try = d_in - widen_step
        s_try = score_bounds(d_in_try, d_out)
        if s_try > base + eps:
            d_in, base = d_in_try, s_try
            improved = True

        # (c) déplacer OUT seul
        d_out_try = d_out + widen_step
        s_try = score_bounds(d_in, d_out_try)
        if s_try > base + eps:
            d_out, base = d_out_try, s_try
            improved = True

    half_opt = 0.5 * (d_out - d_in)
    # score final recompté sur la fenêtre optimale
    score_final = _score_window(residues, t, c_star, half_opt, score_fn)

    return (score_final, n, c_star, half_opt)