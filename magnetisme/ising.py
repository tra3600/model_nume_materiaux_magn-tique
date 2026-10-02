"""Modèle d'Ising 2D : simulation Monte-Carlo (Metropolis séquentiel, damier numpy, Wolff).

Convention (celle du sujet Mines 2022, éq. 3)::

    E = -J * sum_{<ij>} s_i s_j  -  B * sum_i s_i        s_i = ±1

avec conditions aux limites périodiques, kB = 1 et un réseau carré L x L.
Pour J = 1, B = 0 la transition de Curie se produit à Tc = 2 / ln(1 + sqrt 2) ≈ 2.269.

Si numba est installé, les noyaux sont compilés (≈ 50-100x plus rapides) ; sinon un repli
numpy (damier vectorisé) est utilisé pour Metropolis. L'algorithme de Wolff est alors en
Python pur (correct mais lent) : installer numba est recommandé.
"""
from __future__ import annotations

import math

import numpy as np

try:  # numba est optionnel
    from numba import njit
    HAS_NUMBA = True
except ImportError:  # pragma: no cover
    HAS_NUMBA = False

TC_ONSAGER = 2.0 / math.log(1.0 + math.sqrt(2.0))


def _jit(f):
    return njit(cache=True)(f) if HAS_NUMBA else f


# ----------------------------------------------------------------------------
# Noyaux (compilés par numba s'il est présent)
# ----------------------------------------------------------------------------
@_jit
def _seed_numba(seed):
    np.random.seed(seed)


@_jit
def _metropolis_seq(s, T, J, B, nsweeps):
    """Metropolis à sites tirés au hasard (comme monte_carlo du sujet) : nsweeps*L² tests."""
    L = s.shape[0]
    for _ in range(nsweeps * L * L):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)
        nb = (int(s[(i + 1) % L, j]) + int(s[(i - 1) % L, j])
              + int(s[i, (j + 1) % L]) + int(s[i, (j - 1) % L]))
        dE = 2.0 * s[i, j] * (J * nb + B)
        if dE <= 0.0 or np.random.random() < np.exp(-dE / T):
            s[i, j] = -s[i, j]


@_jit
def _wolff_sweep(s, T, J):
    """Retourne des amas de Wolff jusqu'à avoir basculé >= L² spins. Renvoie le nb d'amas."""
    L = s.shape[0]
    p = 1.0 - np.exp(-2.0 * J / T)
    stack = np.empty(L * L, np.int64)
    flipped = 0
    nclusters = 0
    while flipped < L * L:
        i0 = np.random.randint(0, L)
        j0 = np.random.randint(0, L)
        spin = s[i0, j0]
        s[i0, j0] = -spin
        stack[0] = i0 * L + j0
        sp = 1
        flipped += 1
        while sp > 0:
            sp -= 1
            k = stack[sp]
            i = k // L
            j = k % L
            for d in range(4):
                if d == 0:
                    ni = (i + 1) % L
                    nj = j
                elif d == 1:
                    ni = (i - 1) % L
                    nj = j
                elif d == 2:
                    ni = i
                    nj = (j + 1) % L
                else:
                    ni = i
                    nj = (j - 1) % L
                if s[ni, nj] == spin and np.random.random() < p:
                    s[ni, nj] = -spin
                    stack[sp] = ni * L + nj
                    sp += 1
                    flipped += 1
        nclusters += 1
    return nclusters


# ----------------------------------------------------------------------------
# Classe principale
# ----------------------------------------------------------------------------
class Ising2D:
    """Réseau carré L x L de spins ±1 (L pair), conditions périodiques.

    Paramètres
    ----------
    L : côté du réseau (pair, pour le damier et l'ordre antiferro).
    T : température en unités J/kB.
    J : intégrale d'échange (>0 ferro, <0 antiferro).
    B : champ extérieur (couplage -B*sum s_i).
    init : 'froid' (tous +1), 'chaud' (aléatoire), 'neel' (damier) ou tableau LxL.
    seed : graine reproductible.
    """

    def __init__(self, L=32, T=2.0, J=1.0, B=0.0, init="froid", seed=None):
        if L % 2:
            raise ValueError("L doit être pair")
        self.L, self.T, self.J, self.B = int(L), float(T), float(J), float(B)
        self.rng = np.random.default_rng(seed)
        if HAS_NUMBA:
            _seed_numba(int(self.rng.integers(2**31 - 1)))
        if isinstance(init, np.ndarray):
            self.s = init.astype(np.int8).copy()
        elif init == "froid":
            self.s = np.ones((L, L), np.int8)
        elif init == "chaud":
            self.s = self.rng.choice(np.array([-1, 1], np.int8), size=(L, L))
        elif init == "neel":
            self.s = neel(L)
        else:
            raise ValueError(f"init inconnu : {init!r}")
        ii, jj = np.indices((L, L))
        self._parite = (ii + jj) % 2

    # -- dynamique ---------------------------------------------------------
    def balayage(self, n=1, algo="metropolis"):
        """Avance de n balayages (1 balayage = L² tentatives, ou ≥ L² spins pour Wolff)."""
        if algo == "wolff":
            if self.B != 0.0 or self.J <= 0.0:
                raise ValueError("Wolff : requiert B = 0 et J > 0")
            for _ in range(n):
                self._wolff()
            return
        if algo != "metropolis":
            raise ValueError(f"algo inconnu : {algo!r}")
        if HAS_NUMBA:
            _metropolis_seq(self.s, self.T, self.J, self.B, int(n))
        else:
            for _ in range(n):
                self._damier()

    def _damier(self):
        """Un balayage Metropolis vectorisé en damier (repli sans numba)."""
        s = self.s
        for c in (0, 1):
            nb = (np.roll(s, 1, 0) + np.roll(s, -1, 0)
                  + np.roll(s, 1, 1) + np.roll(s, -1, 1)).astype(np.int16)
            dE = 2.0 * s * (self.J * nb + self.B)
            u = self.rng.random(s.shape)
            ok = (dE <= 0) | (u < np.exp(-np.clip(dE, 0, None) / self.T))
            s[(self._parite == c) & ok] *= -1

    def _wolff(self):
        if HAS_NUMBA:
            return _wolff_sweep(self.s, self.T, self.J)
        return _wolff_python(self.s, self.T, self.J, self.rng)  # pragma: no cover

    # -- observables -------------------------------------------------------
    def aimantation(self):
        """m = (1/n) sum s_i  (définition du sujet, question 18)."""
        return float(self.s.mean())

    def aimantation_alternee(self):
        """Paramètre d'ordre antiferro : (1/n) sum (-1)^(i+j) s_i."""
        return float((self.s * (1 - 2 * self._parite)).mean())

    def energie(self):
        """Énergie par spin."""
        s = self.s.astype(np.int16)
        liens = (s * np.roll(s, -1, 0)).sum() + (s * np.roll(s, -1, 1)).sum()
        return float(-self.J * liens / s.size - self.B * s.mean())

    def energie_totale(self):
        return self.energie() * self.L ** 2

    def copie(self):
        c = Ising2D.__new__(Ising2D)
        c.__dict__.update(self.__dict__)
        c.s = self.s.copy()
        return c


def neel(L):
    """Configuration antiferromagnétique en damier (+1/-1 alternés), cf. figure 3 du sujet."""
    ii, jj = np.indices((L, L))
    return (1 - 2 * ((ii + jj) % 2)).astype(np.int8)


def _wolff_python(s, T, J, rng):  # pragma: no cover - repli lent sans numba
    L = s.shape[0]
    p = 1.0 - math.exp(-2.0 * J / T)
    flipped = 0
    n = 0
    while flipped < L * L:
        i0, j0 = int(rng.integers(L)), int(rng.integers(L))
        spin = s[i0, j0]
        s[i0, j0] = -spin
        pile = [(i0, j0)]
        flipped += 1
        while pile:
            i, j = pile.pop()
            for ni, nj in (((i + 1) % L, j), ((i - 1) % L, j), (i, (j + 1) % L), (i, (j - 1) % L)):
                if s[ni, nj] == spin and rng.random() < p:
                    s[ni, nj] = -spin
                    pile.append((ni, nj))
                    flipped += 1
        n += 1
    return n


# ----------------------------------------------------------------------------
# Mesures & statistiques
# ----------------------------------------------------------------------------
def mesurer(modele: Ising2D, n_therm=500, n_mesures=2000, espacement=1, algo="metropolis"):
    """Thermalise puis enregistre (m, e, m_alt) après chaque `espacement` balayages."""
    modele.balayage(n_therm, algo)
    m = np.empty(n_mesures)
    e = np.empty(n_mesures)
    ma = np.empty(n_mesures)
    for k in range(n_mesures):
        modele.balayage(espacement, algo)
        m[k], e[k], ma[k] = modele.aimantation(), modele.energie(), modele.aimantation_alternee()
    return {"m": m, "e": e, "m_alt": ma}


def jackknife(fonction, series, n_blocs=20):
    """Estimation et erreur jackknife de fonction(*series) sur des blocs contigus."""
    series = [np.asarray(x) for x in series]
    n = (len(series[0]) // n_blocs) * n_blocs
    series = [x[:n] for x in series]
    plein = fonction(*series)
    mask = np.arange(n) // (n // n_blocs)
    vals = np.array([fonction(*[x[mask != b] for x in series]) for b in range(n_blocs)])
    err = np.sqrt((n_blocs - 1) * np.mean((vals - vals.mean()) ** 2))
    return float(plein), float(err)


def thermodynamique(mes: dict, L: int, T: float, n_blocs=20, ordre="ferro"):
    """Moyennes <|m|>, <e>, chi, C, Binder U4 avec barres d'erreur jackknife."""
    m = mes["m"] if ordre == "ferro" else mes["m_alt"]
    e = mes["e"]
    n = L * L
    res = {}
    res["M"] = jackknife(lambda a: np.mean(np.abs(a)), [m], n_blocs)
    res["E"] = jackknife(lambda a: np.mean(a), [e], n_blocs)
    res["chi"] = jackknife(lambda a: n * (np.mean(a ** 2) - np.mean(np.abs(a)) ** 2) / T, [m], n_blocs)
    res["C"] = jackknife(lambda a: n * (np.mean(a ** 2) - np.mean(a) ** 2) / T ** 2, [e], n_blocs)
    res["U4"] = jackknife(lambda a: 1 - np.mean(a ** 4) / (3 * np.mean(a ** 2) ** 2), [m], n_blocs)
    return res


def temps_autocorrelation(x, c=6.0):
    """Temps d'autocorrélation intégré (fenêtre automatique de Sokal), en unités d'échantillons."""
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    if acf[0] == 0:
        return 0.5
    acf /= acf[0]
    tau = np.cumsum(acf) - 0.5
    for w in range(1, n):
        if w >= c * tau[w]:
            return float(max(tau[w], 0.5))
    return float(max(tau[-1], 0.5))
