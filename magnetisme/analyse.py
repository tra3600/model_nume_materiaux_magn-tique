"""Expériences numériques : balayages en température, mise à l'échelle, hystérésis, trempe."""
from __future__ import annotations

import numpy as np

from .ising import Ising2D, mesurer, thermodynamique, temps_autocorrelation
from .theorie import TC_ONSAGER


def balayage_temperature(L, temperatures, n_therm=1000, n_mesures=4000, J=1.0, ordre="ferro",
                         algo="metropolis", seed=0):
    """Grandeurs thermodynamiques vs T (départ de l'état ordonné adapté, thermalisation puis mesures).

    Renvoie un dict {grandeur: array (n_T, 2)} avec colonnes (valeur, erreur jackknife).
    """
    init = "froid" if ordre == "ferro" else "neel"
    sortie = {k: [] for k in ("M", "E", "chi", "C", "U4")}
    for k, T in enumerate(temperatures):
        mod = Ising2D(L, T, J=J, init=init, seed=seed + k)
        res = thermodynamique(mesurer(mod, n_therm, n_mesures, algo=algo), L, T, ordre=ordre)
        for g in sortie:
            sortie[g].append(res[g])
    out = {g: np.array(v) for g, v in sortie.items()}
    out["T"] = np.asarray(temperatures, float)
    return out


def ajuster_exposant(x, y, xmin=None, xmax=None):
    """Pente (et ordonnée) d'un ajustement en loi de puissance y = A x^p sur [xmin, xmax]."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = (y > 0) & (x > 0)
    if xmin is not None:
        m &= x >= xmin
    if xmax is not None:
        m &= x <= xmax
    p, lnA = np.polyfit(np.log(x[m]), np.log(y[m]), 1)
    return float(p), float(np.exp(lnA))


def croisement_binder(T, U_petit, U_grand):
    """Température où les cumulants de Binder de deux tailles se croisent (interp. linéaire)."""
    d = np.asarray(U_grand) - np.asarray(U_petit)
    for k in range(len(d) - 1):
        if d[k] == 0:
            return float(T[k])
        if d[k] * d[k + 1] < 0:
            return float(T[k] - d[k] * (T[k + 1] - T[k]) / (d[k + 1] - d[k]))
    return float("nan")


def boucle_hysteresis(L=48, T=1.5, Bmax=1.5, n_champs=61, balayages_par_champ=40, cycles=1, seed=0):
    """Cycle d'hystérésis M(B) : on descend de +Bmax à -Bmax puis on remonte (dynamique Metropolis).

    Renvoie (B, m) concaténés ; la boucle ouverte reflète la métastabilité (nucléation de domaines).
    """
    mod = Ising2D(L, T, B=Bmax, init="froid", seed=seed)
    mod.balayage(100)
    champs = np.concatenate([np.linspace(Bmax, -Bmax, n_champs), np.linspace(-Bmax, Bmax, n_champs)])
    B, m = [], []
    for _ in range(cycles):
        for b in champs:
            mod.B = float(b)
            mod.balayage(balayages_par_champ // 2)
            acc = []
            for _ in range(balayages_par_champ // 2):
                mod.balayage(1)
                acc.append(mod.aimantation())
            B.append(b)
            m.append(np.mean(acc))
    return np.array(B), np.array(m)


def trempe(L=256, T=0.8, temps=(1, 3, 10, 30, 100, 300, 1000), seed=0):
    """Trempe depuis T=∞ (config aléatoire) vers T<Tc : croissance de domaines (coarsening).

    Renvoie ({t: configuration}, t, longueur_caractéristique) avec ℓ(t) = 1/densité de parois,
    attendue ∝ t^(1/2) (Lifshitz-Allen-Cahn).
    """
    mod = Ising2D(L, T, init="chaud", seed=seed)
    instantanes, ts, ell = {}, [], []
    t_prec = 0
    for t in sorted(set(int(x) for x in temps)):
        mod.balayage(t - t_prec)
        t_prec = t
        s = mod.s
        rho = ((s != np.roll(s, 1, 0)).sum() + (s != np.roll(s, 1, 1)).sum()) / (2 * s.size)
        ts.append(t)
        ell.append(1.0 / rho)
        instantanes[t] = s.copy()
    return instantanes, np.array(ts), np.array(ell)


def ralentissement_critique(tailles=(8, 16, 32), n_mesures=6000, seed=0):
    """τ_int de |m| à Tc pour Metropolis et Wolff (en balayages) : exposant dynamique z."""
    res = {"L": np.array(tailles), "metropolis": [], "wolff": []}
    for L in tailles:
        for algo in ("metropolis", "wolff"):
            mod = Ising2D(L, TC_ONSAGER, init="chaud", seed=seed)
            mes = mesurer(mod, 500, n_mesures, algo=algo)
            res[algo].append(temps_autocorrelation(np.abs(mes["m"])))
    res["metropolis"], res["wolff"] = np.array(res["metropolis"]), np.array(res["wolff"])
    return res
