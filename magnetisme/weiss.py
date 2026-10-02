"""Domaines de Weiss : composantes connexes de spins identiques (conditions périodiques)."""
from __future__ import annotations

import numpy as np
from scipy import ndimage


def domaines(s: np.ndarray):
    """Étiquette chaque spin par son domaine de Weiss (4-connexité, bords périodiques).

    Renvoie (etiquettes, n_domaines) avec des numéros 0..n-1 ; version vectorisée (scipy) de
    construire_domaines_weiss du sujet (question 24).
    """
    L = s.shape[0]
    lab = np.zeros(s.shape, np.int64)
    n = 0
    for val in (1, -1):
        l, k = ndimage.label(s == val)
        lab[s == val] = l[s == val] + n
        n += k
    parent = np.arange(n + 1)

    def trouver(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    # fusion des domaines coupés par le bord périodique
    paires = [(lab[0, j], lab[-1, j]) for j in range(L) if s[0, j] == s[-1, j]]
    paires += [(lab[i, 0], lab[i, -1]) for i in range(L) if s[i, 0] == s[i, -1]]
    for a, b in paires:
        ra, rb = trouver(a), trouver(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    racines = np.array([trouver(a) for a in range(n + 1)])
    _, inv = np.unique(racines[lab], return_inverse=True)
    inv = inv.reshape(s.shape)
    # numérotation dans l'ordre de lecture (le spin 0 est dans le domaine 0), comme le sujet
    _, premier = np.unique(inv.ravel(), return_index=True)
    ordre = np.argsort(np.argsort(premier))
    inv = ordre[inv]
    return inv, int(inv.max()) + 1


def tailles_domaines(etiquettes: np.ndarray) -> np.ndarray:
    """Surface (nombre de spins) de chaque domaine."""
    return np.bincount(etiquettes.ravel())


def statistiques(s: np.ndarray) -> dict:
    """Nombre de domaines, plus grand domaine (fraction), taille moyenne, densité de parois."""
    lab, n = domaines(s)
    tailles = tailles_domaines(lab)
    N = s.size
    parois = (s != np.roll(s, 1, 0)).sum() + (s != np.roll(s, 1, 1)).sum()
    return {
        "n_domaines": n,
        "plus_grand": tailles.max() / N,
        "taille_moyenne": N / n,
        "densite_parois": parois / (2 * N),
    }


def correlation(s: np.ndarray, rmax=None):
    """Fonction de corrélation connexe G(r) = <s_0 s_r> - <s>² (moyenne sur les deux axes, FFT)."""
    L = s.shape[0]
    x = s.astype(float)
    f = np.fft.fft2(x)
    c = np.fft.ifft2(f * np.conj(f)).real / x.size
    c -= x.mean() ** 2
    r = np.arange(L // 2 + 1)
    g = 0.5 * (c[r, 0] + c[0, r])
    return (r, g) if rmax is None else (r[:rmax], g[:rmax])


def longueur_correlation(r, g):
    """ξ par ajustement exponentiel G(r) ≈ A exp(-r/ξ) sur la zone où G>0 et r>=1."""
    mask = (r >= 1) & (g > 1e-3 * g[0])
    if mask.sum() < 3:
        return float("nan")
    pente = np.polyfit(r[mask], np.log(g[mask]), 1)[0]
    return float(-1 / pente) if pente < 0 else float("nan")
