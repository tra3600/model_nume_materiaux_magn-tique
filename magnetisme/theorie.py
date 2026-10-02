"""Théories analytiques : champ moyen de Weiss, Brillouin, Landau, solution exacte d'Onsager.

Unités : kB = 1, J = 1 sauf mention contraire.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import optimize, special

TC_ONSAGER = 2.0 / math.log(1.0 + math.sqrt(2.0))      # ≈ 2.2692
TC_CHAMP_MOYEN_2D = 4.0                                 # z*J avec z = 4 voisins

# Exposants critiques : (champ moyen, Ising 2D exact)
EXPOSANTS = {
    "beta":  (1 / 2, 1 / 8),    # M ~ (Tc-T)^beta
    "gamma": (1.0, 7 / 4),      # chi ~ |T-Tc|^-gamma
    "nu":    (1 / 2, 1.0),      # xi ~ |T-Tc|^-nu
    "alpha": (0.0, 0.0),        # C ~ saut (MF) / log (2D)
    "delta": (3.0, 15.0),       # M ~ B^(1/delta) à Tc
}


# ---------------------------------------------------------------- champ moyen (sujet partie I)
def f(x, t):
    """Question 2 : f(x,t) = tanh(x/t) - x ; son zéro est l'aimantation réduite m = tanh(m/t)."""
    return np.tanh(x / t) - x


def dicho(fonction, t, a, b, eps):
    """Question 3 : zéro de fonction(., t) sur [a,b] à eps près (dichotomie, O(log2((b-a)/eps)))."""
    fa = fonction(a, t)
    while (b - a) / 2.0 > eps:
        m = (a + b) / 2.0
        fm = fonction(m, t)
        if fm == 0:
            return m
        if fa * fm < 0:
            b = m
        else:
            a, fa = m, fm
    return (a + b) / 2.0


def aimantation_reduite(t, b=0.0):
    """m(t, b) de l'équation d'état m = tanh((m+b)/t) (t = T/Tc, b = B/(zJ)).

    Pour b = 0 : branche stable m>0 (0 si t >= 1). Pour b != 0 : racine de même signe que b
    (la plus stable), trouvée par brentq.
    """
    t = float(t)
    if b == 0.0:
        if t >= 1.0:
            return 0.0
        return float(optimize.brentq(lambda x: f(x, t), 1e-9 if t > 0.05 else 1e-300, 1.0, xtol=1e-14))
    g = lambda x: math.tanh((x + b) / t) - x
    return float(optimize.brentq(g, 0.0, 1.0) if b > 0 else optimize.brentq(g, -1.0, 0.0))


def construction_liste_m(t1, t2, n=500, eps=1e-6):
    """Question 5 : n solutions m(t) pour t linéaire de t1 à t2 (bornes incluses)."""
    ts = np.linspace(t1, t2, n)
    return [0.0 if t >= 1 else dicho(f, t, 0.001, 1.0, eps) for t in ts]


def susceptibilite_champ_moyen(t):
    """Loi de Curie-Weiss réduite : chi = 1/(t-1) pour t>1 ; pour t<1 : (1-m²)/(t-(1-m²))."""
    t = np.asarray(t, float)
    out = np.empty_like(t)
    for k, tk in np.ndenumerate(t):
        if tk == 1:
            out[k] = np.inf
        elif tk > 1:
            out[k] = 1.0 / (tk - 1.0)
        else:
            m = aimantation_reduite(tk)
            q = 1 - m * m
            out[k] = q / (tk - q)
    return out


def landau(m, t):
    """Énergie libre de champ moyen réduite f(m) = m²/2 - t ln(2 cosh(m/t))  (minimum en m=tanh(m/t))."""
    m = np.asarray(m, float)
    return m ** 2 / 2 - t * (np.log(2) + np.abs(m / t) + np.log1p(np.exp(-2 * np.abs(m / t))) - np.abs(m / t))


def brillouin(x, J):
    """Fonction de Brillouin B_J(x) (J = moment cinétique total, demi-entier)."""
    x = np.asarray(x, float)
    a, b = (2 * J + 1) / (2 * J), 1 / (2 * J)
    xs = np.where(np.abs(x) < 1e-8, 1e-8, x)
    out = a / np.tanh(a * xs) - b / np.tanh(b * xs)
    return np.where(np.abs(x) < 1e-8, (J + 1) / (3 * J) * x, out)


def aimantation_brillouin(t, J):
    """M/Ms(T/Tc) de la théorie de Weiss quantique : m = B_J( 3J/(J+1) m/t ). J=∞ → Langevin, J=1/2 → tanh."""
    t = np.atleast_1d(np.asarray(t, float))
    out = np.zeros_like(t)
    c = 3.0 * J / (J + 1.0)
    for k, tk in enumerate(t):
        if tk >= 1:
            continue
        g = lambda m: float(brillouin(c * m / tk, J)) - m
        out[k] = optimize.brentq(g, 1e-6, 1.0, xtol=1e-12) if g(1e-6) > 0 else 0.0
    return out


# ---------------------------------------------------------------- Onsager (Ising 2D exact, B=0)
def aimantation_onsager(T, J=1.0):
    """Aimantation spontanée exacte (Yang 1952) : M = (1 - sinh(2J/T)^-4)^(1/8) pour T<Tc."""
    T = np.asarray(T, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = 1.0 - np.sinh(2 * J / T) ** -4.0
    return np.where(T < TC_ONSAGER * abs(J), np.clip(x, 0, None) ** 0.125, 0.0)


def _k(T, J):
    K = J / T
    return 2 * np.sinh(2 * K) / np.cosh(2 * K) ** 2


def energie_onsager(T, J=1.0):
    """Énergie interne par spin exacte u(T) = -J coth(2K)[1 + 2/π (2 tanh²2K - 1) K1(k)]."""
    T = np.asarray(T, float)
    K = J / T
    k = _k(T, J)
    th = np.tanh(2 * K)
    return -J / th * (1 + (2 / np.pi) * (2 * th ** 2 - 1) * special.ellipk(k ** 2))


def chaleur_specifique_onsager(T, J=1.0):
    """C/kB par spin exacte (divergence logarithmique en Tc)."""
    T = np.asarray(T, float)
    K = J / T
    k = _k(T, J)
    th = np.tanh(2 * K)
    K1, E1 = special.ellipk(k ** 2), special.ellipe(k ** 2)
    return (1 / np.pi) * (2 * K / th) ** 2 * (K1 - E1 - (1 - th ** 2) * (np.pi / 2 + (2 * th ** 2 - 1) * K1))


def longueur_correlation_haute_temperature(T, J=1.0):
    """ξ(T>Tc) exacte de la fonction de corrélation le long d'un axe : ξ^-1 = 2K + ln tanh K."""
    K = J / np.asarray(T, float)
    return 1.0 / (2 * K + np.log(np.tanh(K)))
