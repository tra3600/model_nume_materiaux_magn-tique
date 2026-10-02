"""Génération de toutes les figures (PNG + GIF) dans figures/.

`generer_tout(rapide=True)` utilise des tailles/durées réduites (tests, démonstration < 1 min).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use(matplotlib.get_backend() if matplotlib.get_backend() != "agg" else "Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

from . import analyse, materiaux, theorie, weiss  # noqa: E402
from .ising import Ising2D, mesurer  # noqa: E402
from .theorie import TC_ONSAGER  # noqa: E402

# Palette daltonien-compatible (Okabe-Ito)
C = dict(bleu="#0072B2", orange="#E69F00", vert="#009E73", rouge="#D55E00",
         violet="#CC79A7", ciel="#56B4E9", noir="#222222", gris="#888888")
SPINS = ListedColormap(["#f2f2f2", "#1b2a49"])      # -1 clair, +1 foncé (cf. figure 5 du sujet)
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 150, "savefig.bbox": "tight", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.alpha": 0.25, "axes.prop_cycle": plt.cycler(color=list(C.values())[:7]),
})


def _sauver(fig, dossier: Path, nom: str) -> Path:
    dossier.mkdir(parents=True, exist_ok=True)
    chemin = dossier / nom
    fig.savefig(chemin)
    plt.close(fig)
    return chemin


def _params(rapide):
    if rapide:
        return dict(L=16, tailles=(8, 16), nT=15, therm=200, mes=600, L_snap=48, L_trempe=64,
                    L_hys=24, n_champs=21, t_max=300)
    return dict(L=32, tailles=(8, 16, 32, 64), nT=41, therm=2000, mes=10000, L_snap=100, L_trempe=256,
                L_hys=64, n_champs=61, t_max=3000)


# ------------------------------------------------------------------------------------------------
def fig_champ_moyen(dossier, rapide=False):
    """Fig. 1 du sujet (m(t)) comparée à la solution exacte 2D et à la simulation."""
    p = _params(rapide)
    t = np.linspace(0.01, 2.0, 400)
    m_mf = np.array([theorie.aimantation_reduite(x) for x in t])
    T = t * theorie.TC_CHAMP_MOYEN_2D
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(t, m_mf, color=C["bleu"], lw=2.2, label="Champ moyen  $m=\\tanh(m/t)$")
    ax[0].plot(t, theorie.aimantation_onsager(t * TC_ONSAGER), color=C["rouge"], lw=2.2,
               label="Ising 2D exact (Onsager–Yang)")
    ts = np.linspace(0.3, 1.6, 14 if not rapide else 8)
    r = analyse.balayage_temperature(p["L"] * 2, ts * TC_ONSAGER, p["therm"], p["mes"])
    ax[0].errorbar(ts, r["M"][:, 0], r["M"][:, 1], fmt="o", ms=4, color=C["noir"],
                   label=f"Monte-Carlo ($L={p['L']*2}$)")
    ax[0].axvline(1, color=C["gris"], ls=":")
    ax[0].set(xlabel="$t=T/T_c$", ylabel="aimantation réduite $m$",
              title="Transition ferro → paramagnétique")
    ax[0].legend(frameon=False)
    ax[0].annotate("ferro", (0.35, 0.15), color=C["gris"])
    ax[0].annotate("para", (1.45, 0.15), color=C["gris"])
    # exposants : log-log près de Tc
    eps = np.logspace(-3, -0.7, 60)
    ax[1].loglog(eps, [theorie.aimantation_reduite(1 - e) for e in eps], color=C["bleu"], lw=2,
                 label="champ moyen, $\\beta=1/2$")
    ax[1].loglog(eps, theorie.aimantation_onsager(TC_ONSAGER * (1 - eps)), color=C["rouge"], lw=2,
                 label="Ising 2D, $\\beta=1/8$")
    ax[1].set(xlabel="$1-T/T_c$", ylabel="$m$", title="Exposant critique $\\beta$ : $m\\propto(1-t)^{\\beta}$")
    ax[1].legend(frameon=False)
    return _sauver(fig, dossier, "01_aimantation_champ_moyen.png")


def fig_instantanes(dossier, rapide=False):
    """Fig. 5 du sujet : configurations à différentes températures."""
    p = _params(rapide)
    Ts = [1.0, 2.0, TC_ONSAGER, 2.6, 4.0]
    fig, ax = plt.subplots(1, len(Ts), figsize=(3.1 * len(Ts), 3.5))
    for a, T in zip(ax, Ts):
        mod = Ising2D(p["L_snap"], T, init="chaud", seed=3)
        mod.balayage(300 if rapide else 1500, algo="wolff" if T > 2 else "metropolis")
        a.imshow(mod.s, cmap=SPINS, vmin=-1, vmax=1, interpolation="nearest")
        a.set_title(f"T = {T:.3f}" + ("  ($T_c$)" if abs(T - TC_ONSAGER) < 1e-9 else "") +
                    f"\n$|m|$ = {abs(mod.aimantation()):.2f}")
        a.set_xticks([]); a.set_yticks([]); a.grid(False)
    fig.suptitle("Évolution des domaines avec la température (noir : spin ↑, clair : spin ↓)", y=1.02)
    return _sauver(fig, dossier, "02_instantanes_temperature.png")


def fig_observables(dossier, rapide=False):
    p = _params(rapide)
    T = np.linspace(1.2, 3.6, p["nT"])
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    Tf = np.linspace(1.2, 3.6, 400)
    for k, L in enumerate(p["tailles"]):
        r = analyse.balayage_temperature(L, T, p["therm"], p["mes"], seed=100 * k)
        col = list(C.values())[k]
        for a, g in zip(ax.ravel(), ("M", "E", "chi", "C")):
            a.errorbar(T, r[g][:, 0], r[g][:, 1], fmt="o-", ms=3, lw=1, color=col, label=f"$L={L}$")
    ax[0, 0].plot(Tf, theorie.aimantation_onsager(Tf), "k--", lw=1.5, label="exact ($L=\\infty$)")
    ax[0, 1].plot(Tf, theorie.energie_onsager(Tf), "k--", lw=1.5)
    Cf = theorie.chaleur_specifique_onsager(Tf[np.abs(Tf - TC_ONSAGER) > 0.01])
    ax[1, 1].plot(Tf[np.abs(Tf - TC_ONSAGER) > 0.01], Cf, "k--", lw=1.5)
    for a, yl, ti in zip(ax.ravel(), ("$\\langle|m|\\rangle$", "énergie par spin $e$",
                                       "susceptibilité $\\chi$", "chaleur spécifique $C$"),
                         ("Aimantation", "Énergie", "Susceptibilité", "Chaleur spécifique")):
        a.axvline(TC_ONSAGER, color=C["gris"], ls=":")
        a.set(xlabel="$T$ (en $J/k_B$)", ylabel=yl, title=ti)
    ax[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Ising 2D : observables Monte-Carlo vs solution exacte (pointillés : $T_c=2.269$)")
    fig.tight_layout()
    return _sauver(fig, dossier, "03_observables.png")


def fig_binder_fss(dossier, rapide=False):
    p = _params(rapide)
    T = np.linspace(2.0, 2.55, 23 if not rapide else 11)
    res = {}
    for k, L in enumerate(p["tailles"]):
        res[L] = analyse.balayage_temperature(L, T, p["therm"], p["mes"], seed=7 * k + 1)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for k, L in enumerate(p["tailles"]):
        col = list(C.values())[k]
        ax[0].errorbar(T, res[L]["U4"][:, 0], res[L]["U4"][:, 1], fmt="o-", ms=3, color=col, label=f"$L={L}$")
        x = (T - TC_ONSAGER) * L
        ax[1].plot(x, res[L]["M"][:, 0] * L ** 0.125, "o", ms=4, color=col, label=f"$L={L}$")
        ax[2].plot(x, res[L]["chi"][:, 0] / L ** 1.75, "o", ms=4, color=col)
    Ls = list(p["tailles"])
    Tx = analyse.croisement_binder(T, res[Ls[-2]]["U4"][:, 0], res[Ls[-1]]["U4"][:, 0])
    ax[0].axvline(TC_ONSAGER, color=C["gris"], ls=":")
    ax[0].axhline(0.6107, color=C["gris"], ls=":")
    ax[0].set(xlabel="$T$", ylabel="$U_4=1-\\langle m^4\\rangle/3\\langle m^2\\rangle^2$",
              title=f"Cumulant de Binder — croisement $L={Ls[-2]},{Ls[-1]}$ : $T_x$ = {Tx:.3f}\n"
                    f"(exact $T_c$ = {TC_ONSAGER:.3f},  $U_4^*\\approx0.611$)")
    ax[0].legend(frameon=False)
    ax[1].set(xlabel="$(T-T_c)\\,L^{1/\\nu}$", ylabel="$\\langle|m|\\rangle L^{\\beta/\\nu}$",
              title="Collapse : $\\beta=1/8,\\ \\nu=1$")
    ax[2].set(xlabel="$(T-T_c)\\,L^{1/\\nu}$", ylabel="$\\chi\\,L^{-\\gamma/\\nu}$",
              title="Collapse : $\\gamma=7/4,\\ \\nu=1$")
    return _sauver(fig, dossier, "04_binder_finite_size_scaling.png")


def fig_domaines(dossier, rapide=False):
    p = _params(rapide)
    Ts = [1.8, 2.15, TC_ONSAGER, 2.5]
    fig, ax = plt.subplots(2, len(Ts), figsize=(3.3 * len(Ts), 7))
    rng = np.random.default_rng(0)
    gris = plt.get_cmap("nipy_spectral")
    for k, T in enumerate(Ts):
        mod = Ising2D(p["L_snap"], T, init="chaud", seed=11)
        mod.balayage(400 if rapide else 2500, algo="wolff")
        lab, n = weiss.domaines(mod.s)
        perm = rng.permutation(n)
        ax[0, k].imshow(gris(perm[lab] / max(n - 1, 1)), interpolation="nearest")
        st = weiss.statistiques(mod.s)
        ax[0, k].set_title(f"T = {T:.3f}\n{st['n_domaines']} domaines, plus grand : {st['plus_grand']:.0%}")
        ax[0, k].set_xticks([]); ax[0, k].set_yticks([]); ax[0, k].grid(False)
        tailles = weiss.tailles_domaines(lab)
        tailles = tailles[tailles > 0]
        b = np.unique(np.logspace(0, np.log10(tailles.max() + 1), 15).astype(int))
        h, e = np.histogram(tailles, bins=b)
        ax[1, k].loglog(np.sqrt(e[:-1] * e[1:]), h / np.diff(e), "o-", color=C["bleu"])
        ax[1, k].set(xlabel="taille du domaine (spins)", ylabel="densité $P(s)$")
    ax[1, 0].set_title("distribution des tailles", fontsize=9)
    fig.suptitle("Domaines de Weiss (composantes connexes, fig. 6-7 du sujet) — à $T_c$ : loi de puissance", y=1.0)
    fig.tight_layout()
    return _sauver(fig, dossier, "05_domaines_weiss.png")


def fig_correlations(dossier, rapide=False):
    p = _params(rapide)
    Ts = [2.4, 2.6, 3.0, 3.6]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    xis = []
    for T in Ts:
        L = p["L_snap"]
        Gs = []
        mod = Ising2D(L, T, init="chaud", seed=5)
        mod.balayage(500)
        for _ in range(60 if not rapide else 10):
            mod.balayage(5)
            r, g = weiss.correlation(mod.s)
            Gs.append(g)
        g = np.mean(Gs, 0)
        ax[0].semilogy(r, np.clip(g, 1e-4, None), "o-", ms=3, label=f"T={T}")
        xis.append(weiss.longueur_correlation(r, g))
    ax[0].set(xlabel="distance $r$", ylabel="$G(r)=\\langle s_0s_r\\rangle_c$", title="Corrélations (T > Tc)")
    ax[0].legend(frameon=False)
    Tf = np.linspace(2.35, 3.8, 100)
    ax[1].plot(Tf, theorie.longueur_correlation_haute_temperature(Tf), color=C["rouge"], label="exact $\\xi(T)$")
    ax[1].plot(Ts, xis, "o", color=C["noir"], label="Monte-Carlo")
    ax[1].set(xlabel="$T$", ylabel="longueur de corrélation $\\xi$", ylim=(0, 8),
              title="Divergence de $\\xi$ en $T_c$ ($\\nu=1$)")
    ax[1].legend(frameon=False)
    return _sauver(fig, dossier, "06_correlations.png")


def fig_hysteresis(dossier, rapide=False):
    p = _params(rapide)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for k, T in enumerate((1.2, 1.8, 2.2, 2.6)):
        B, m = analyse.boucle_hysteresis(p["L_hys"], T, 1.5, p["n_champs"], 20 if rapide else 60, seed=k)
        ax[0].plot(B, m, color=list(C.values())[k], label=f"T={T}")
    ax[0].set(xlabel="champ $B$", ylabel="aimantation $m$", title="Cycles d'hystérésis (Metropolis)")
    ax[0].legend(frameon=False)
    # instantanés au champ coercitif
    mod = Ising2D(p["L_hys"], 1.5, B=1.0, init="froid", seed=2)
    mod.balayage(50)
    mod.B = -0.3
    mod.balayage(60 if not rapide else 20)
    ax[1].imshow(mod.s, cmap=SPINS, vmin=-1, vmax=1, interpolation="nearest")
    ax[1].set_title("Nucléation de domaines inverses ($T=1.5$, $B=-0.3$)")
    ax[1].set_xticks([]); ax[1].set_yticks([]); ax[1].grid(False)
    return _sauver(fig, dossier, "07_hysteresis.png")


def fig_trempe(dossier, rapide=False):
    p = _params(rapide)
    temps = [t for t in (1, 3, 10, 30, 100, 300, 1000, 3000) if t <= p["t_max"]]
    inst, ts, ell = analyse.trempe(p["L_trempe"], 0.8, temps, seed=1)
    montrer = [t for t in (3, 30, 300, temps[-1]) if t in inst][:4]
    fig = plt.figure(figsize=(14, 4))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 1, 1.6])
    for k, t in enumerate(montrer):
        a = fig.add_subplot(gs[k])
        a.imshow(inst[t], cmap=SPINS, vmin=-1, vmax=1, interpolation="nearest")
        a.set_title(f"t = {t} balayages"); a.set_xticks([]); a.set_yticks([]); a.grid(False)
    a = fig.add_subplot(gs[4])
    a.loglog(ts, ell, "o", color=C["bleu"], label="Monte-Carlo")
    sel = ts >= 10
    p_, A = analyse.ajuster_exposant(ts[sel], ell[sel])
    a.loglog(ts, A * ts ** p_, "-", color=C["rouge"], label=f"ajustement $t^{{{p_:.2f}}}$")
    a.loglog(ts, ell[sel][0] * (ts / ts[sel][0]) ** 0.5, "--", color=C["gris"], label="loi $t^{1/2}$")
    a.set(xlabel="temps $t$ (balayages)", ylabel="taille de domaine $\\ell=1/\\rho_{parois}$",
          title="Coarsening après trempe $T=0.8$")
    a.legend(frameon=False)
    return _sauver(fig, dossier, "08_trempe_coarsening.png")


def fig_antiferro(dossier, rapide=False):
    p = _params(rapide)
    T = np.linspace(0.8, 3.6, p["nT"])
    L = p["L"]
    ferro = analyse.balayage_temperature(L, T, p["therm"], p["mes"], J=+1, ordre="ferro", seed=1)
    anti = analyse.balayage_temperature(L, T, p["therm"], p["mes"], J=-1, ordre="neel", seed=1)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].errorbar(T, ferro["M"][:, 0], ferro["M"][:, 1], fmt="o", ms=3, color=C["bleu"], label="ferro $J=+1$ : $|m|$")
    ax[0].errorbar(T, anti["M"][:, 0], anti["M"][:, 1], fmt="s", ms=3, color=C["rouge"],
                   label="antiferro $J=-1$ : $|m_{alt}|$")
    ax[0].plot(T, theorie.aimantation_onsager(T), "k--", lw=1, label="Onsager")
    ax[0].axvline(TC_ONSAGER, color=C["gris"], ls=":")
    ax[0].set(xlabel="$T$", ylabel="paramètre d'ordre", title="Même transition ($T_N=T_C$, réseau bipartite)")
    ax[0].legend(frameon=False, fontsize=8)
    # susceptibilité uniforme de l'antiferro (m uniforme) : maximum cusp en TN
    res_u = []
    for k, t in enumerate(T):
        mod = Ising2D(L, t, J=-1, init="neel", seed=k)
        mes = mesurer(mod, p["therm"], p["mes"])
        res_u.append(L * L * np.var(mes["m"]) / t)
    ax[1].plot(T, res_u, "o-", ms=3, color=C["vert"])
    ax[1].axvline(TC_ONSAGER, color=C["gris"], ls=":")
    ax[1].set(xlabel="$T$", ylabel="$\\chi_{unif}$", title="Susceptibilité uniforme : pic en $T_N$")
    mod = Ising2D(L, 1.0, J=-1, init="neel", seed=1); mod.balayage(200)
    ax[2].imshow(mod.s, cmap=SPINS, vmin=-1, vmax=1, interpolation="nearest")
    ax[2].set_title("Ordre de Néel à $T=1$ (fig. 3 du sujet)")
    ax[2].set_xticks([]); ax[2].set_yticks([]); ax[2].grid(False)
    return _sauver(fig, dossier, "09_antiferromagnetisme.png")


def fig_curie_weiss(dossier, rapide=False):
    p = _params(rapide)
    T = np.linspace(2.7, 6.0, 12 if not rapide else 6)
    r = analyse.balayage_temperature(p["tailles"][-1] * (2 if not rapide else 1), T, p["therm"], p["mes"])
    chi = r["chi"][:, 0]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(T, 1 / chi, "o", color=C["noir"], label="Monte-Carlo")
    tt = np.linspace(2.0, 6.5, 100)
    ax[0].plot(tt, (tt - 4.0) * 1.0, color=C["bleu"], label="Curie-Weiss (champ moyen) $\\theta=4$")
    ax[0].set(xlabel="$T$", ylabel="$1/\\chi$", ylim=(0, 3),
              title="Loi de Curie–Weiss : $\\chi=C/(T-\\theta)$ ?")
    ax[0].legend(frameon=False, fontsize=8)
    g, A = analyse.ajuster_exposant(T - TC_ONSAGER, chi, 0.4, 1.5)
    ax[1].loglog(T - TC_ONSAGER, chi, "o", color=C["noir"], label="Monte-Carlo")
    ax[1].loglog(T - TC_ONSAGER, A * (T - TC_ONSAGER) ** g, color=C["rouge"], label=f"pente ajustée : $-\\gamma$={g:.2f}")
    ax[1].loglog(T - TC_ONSAGER, chi[0] * ((T - TC_ONSAGER) / (T[0] - TC_ONSAGER)) ** -1.75, "--",
                 color=C["gris"], label="exact $\\gamma=7/4$")
    ax[1].set(xlabel="$T-T_c$", ylabel="$\\chi$", title="Susceptibilité au-dessus de $T_c$")
    ax[1].legend(frameon=False, fontsize=8)
    t = np.linspace(0.2, 2.5, 300)
    ax[2].plot(t, np.minimum(theorie.susceptibilite_champ_moyen(t), 30), color=C["bleu"])
    ax[2].set(xlabel="$t=T/T_c$", ylabel="$\\chi_{MF}$ (réduite)", ylim=(0, 30),
              title="Champ moyen : divergence $1/|t-1|$")
    return _sauver(fig, dossier, "10_curie_weiss.png")


def fig_landau(dossier, rapide=False):
    m = np.linspace(-1.2, 1.2, 400)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for k, t in enumerate((0.5, 0.8, 1.0, 1.3)):
        ax[0].plot(m, theorie.landau(m, t) - theorie.landau(0, t), color=list(C.values())[k], label=f"t={t}")
    ax[0].set(xlabel="aimantation $m$", ylabel="$f(m)-f(0)$", ylim=(-0.2, 0.3),
              title="Paysage d'énergie libre (Landau) : double puits → puits unique")
    ax[0].legend(frameon=False)
    t = np.linspace(0.05, 1.5, 200)
    b = [0.0, 0.02, 0.1]
    for bb, col in zip(b, (C["noir"], C["orange"], C["vert"])):
        mm = [theorie.aimantation_reduite(x, bb) if bb else theorie.aimantation_reduite(x) for x in t]
        ax[1].plot(t, mm, color=col, label=f"$B/zJ$={bb}")
    ax[1].set(xlabel="$t=T/T_c$", ylabel="$m$", title="Effet d'un champ : la transition devient un crossover")
    ax[1].legend(frameon=False)
    return _sauver(fig, dossier, "11_landau_champ.png")


def fig_ralentissement(dossier, rapide=False):
    tailles = (8, 16) if rapide else (8, 16, 32, 64)
    r = analyse.ralentissement_critique(tailles, 1500 if rapide else 8000)
    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    ax.loglog(r["L"], r["metropolis"], "o-", color=C["rouge"], label="Metropolis")
    ax.loglog(r["L"], r["wolff"], "s-", color=C["bleu"], label="Wolff (amas)")
    if len(tailles) > 2:
        z, _ = analyse.ajuster_exposant(r["L"], r["metropolis"])
        ax.set_title(f"Ralentissement critique à $T_c$ : Metropolis $\\tau\\sim L^{{{z:.1f}}}$")
    else:
        ax.set_title("Ralentissement critique à $T_c$")
    ax.set(xlabel="taille $L$", ylabel="$\\tau_{int}$ (balayages)")
    ax.legend(frameon=False)
    return _sauver(fig, dossier, "12_ralentissement_critique_wolff.png")


def fig_materiaux(dossier, rapide=False):
    mats = materiaux.table_materiaux()
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    noms = [m["nom"] for m in mats]
    tc = [m["t_curie"] for m in mats]
    cols = [C["bleu"] if "ferro" == m["ordre"][:5] else C["orange"] if "ferri" in m["ordre"] else C["rouge"]
            for m in mats]
    ax[0].barh(noms[::-1], tc[::-1], color=cols[::-1])
    ax[0].axvline(300, color=C["gris"], ls=":"); ax[0].text(310, 0.2, "ambiante", color=C["gris"], fontsize=8)
    ax[0].set(xlabel="$T_C$ ou $T_N$ (K)", title="Base de données : températures de Curie / Néel")
    t = np.linspace(0.001, 1.2, 300)
    for J, nom, col in ((0.5, "J=1/2 (Fe, Co, Ni)", C["bleu"]), (2.5, "J=5/2 (YIG)", C["orange"]),
                        (3.5, "J=7/2 (Gd, EuO)", C["vert"]), (50, "J→∞ (classique)", C["violet"])):
        ax[1].plot(t, theorie.aimantation_brillouin(t, J), color=col, label=nom)
    ax[1].set(xlabel="$T/T_C$", ylabel="$M/M_s$", title="Weiss–Brillouin : $M/M_s=B_J(\\frac{3J}{J+1}\\,m/t)$")
    ax[1].legend(frameon=False, fontsize=8)
    # courbes M(T) en kelvin pour quelques matériaux réels
    T = np.linspace(1, 1500, 600)
    for m, col in zip([x for x in mats if x["nom"] in ("fer", "cobalt", "nickel", "gadolinium", "YIG")],
                      (C["noir"], C["rouge"], C["vert"], C["violet"], C["orange"])):
        ax[2].plot(T, m["ms_kAm"] * theorie.aimantation_brillouin(T / m["t_curie"], m["j_eff"]), color=col,
                   label=m["nom"])
    ax[2].set(xlabel="$T$ (K)", ylabel="$M_s(T)$ (kA/m, modèle)", title="Aimantation spontanée (champ moyen quantique)")
    ax[2].legend(frameon=False, fontsize=8)
    return _sauver(fig, dossier, "13_materiaux_reels.png")


def gif_trempe(dossier, rapide=False):
    """Animation du coarsening (GIF)."""
    from matplotlib.animation import FuncAnimation, PillowWriter
    L = 64 if rapide else 128
    mod = Ising2D(L, 0.8, init="chaud", seed=4)
    fig, ax = plt.subplots(figsize=(4, 4.3))
    im = ax.imshow(mod.s, cmap=SPINS, vmin=-1, vmax=1, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    titre = ax.set_title("t = 0")
    n_img = 20 if rapide else 80
    t_tot = [0]

    def maj(k):
        pas = 1 + k // 4
        mod.balayage(pas)
        t_tot[0] += pas
        im.set_data(mod.s)
        titre.set_text(f"trempe T=0.8 — t = {t_tot[0]} balayages")
        return im, titre

    ani = FuncAnimation(fig, maj, frames=n_img, blit=False)
    dossier.mkdir(parents=True, exist_ok=True)
    chemin = dossier / "14_trempe_animation.gif"
    ani.save(chemin, writer=PillowWriter(fps=12))
    plt.close(fig)
    return chemin


TOUTES = {
    "champ_moyen": fig_champ_moyen, "instantanes": fig_instantanes, "observables": fig_observables,
    "binder": fig_binder_fss, "domaines": fig_domaines, "correlations": fig_correlations,
    "hysteresis": fig_hysteresis, "trempe": fig_trempe, "antiferro": fig_antiferro,
    "curie_weiss": fig_curie_weiss, "landau": fig_landau, "wolff": fig_ralentissement,
    "materiaux": fig_materiaux, "gif": gif_trempe,
}


def generer_tout(dossier="figures", rapide=False, noms=None, verbeux=True):
    import time
    dossier = Path(dossier)
    sortie = []
    for nom, fn in TOUTES.items():
        if noms and nom not in noms:
            continue
        t0 = time.time()
        chemin = fn(dossier, rapide)
        sortie.append(chemin)
        if verbeux:
            print(f"  {chemin}  ({time.time() - t0:.1f}s)")
    return sortie
