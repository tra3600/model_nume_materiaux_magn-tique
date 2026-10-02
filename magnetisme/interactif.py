"""Simulateur interactif (matplotlib) : curseurs T et B, boutons, courbes M(t) en direct.

Lancer :  python -m magnetisme interactif
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, RadioButtons, Slider

from .figures import SPINS
from .ising import Ising2D
from .theorie import TC_ONSAGER, aimantation_onsager


def construire(L=96, T0=2.0, B0=0.0, seed=0):
    """Construit la figure interactive ; renvoie (fig, animation, modele) sans appeler plt.show()."""
    mod = Ising2D(L, T0, B=B0, init="chaud", seed=seed)
    etat = {"algo": "metropolis", "pause": False, "hist": []}
    fig = plt.figure(figsize=(11, 6.2))
    ax_s = fig.add_axes([0.04, 0.28, 0.42, 0.57])
    ax_m = fig.add_axes([0.56, 0.52, 0.40, 0.33])
    im = ax_s.imshow(mod.s, cmap=SPINS, vmin=-1, vmax=1, interpolation="nearest")
    ax_s.set_xticks([]); ax_s.set_yticks([])
    titre = ax_s.set_title("")
    Tf = np.linspace(0.3, 4.0, 200)
    ax_m.plot(Tf, aimantation_onsager(Tf), "k--", lw=1, label="Onsager (B=0)")
    ax_m.axvline(TC_ONSAGER, color="grey", ls=":")
    pt, = ax_m.plot([], [], "o", color="#D55E00", ms=4, alpha=.5)
    cur, = ax_m.plot([], [], "o", color="#0072B2", ms=9)
    ax_m.set(xlim=(0, 4), ylim=(-1.05, 1.05), xlabel="T", ylabel="m", title="m en fonction de T (trace)")
    sT = Slider(fig.add_axes([0.12, 0.17, 0.35, 0.03]), "T", 0.2, 5.0, valinit=T0)
    sB = Slider(fig.add_axes([0.12, 0.11, 0.35, 0.03]), "B", -1.0, 1.0, valinit=B0)
    sV = Slider(fig.add_axes([0.12, 0.05, 0.35, 0.03]), "vitesse", 1, 20, valinit=2, valstep=1)
    radio = RadioButtons(fig.add_axes([0.56, 0.08, 0.15, 0.28]), ("metropolis", "wolff"))
    b_alea = Button(fig.add_axes([0.75, 0.30, 0.2, 0.05]), "T=∞ (aléatoire)")
    b_ordre = Button(fig.add_axes([0.75, 0.23, 0.2, 0.05]), "Tout ↑ (ordonné)")
    b_neel = Button(fig.add_axes([0.75, 0.16, 0.2, 0.05]), "Damier (Néel)")
    b_pause = Button(fig.add_axes([0.75, 0.09, 0.2, 0.05]), "Pause / Reprendre")
    b_J = Button(fig.add_axes([0.75, 0.02, 0.2, 0.05]), "Ferro ⇄ Antiferro (J→−J)")

    def reinit(mode):
        mod.s = Ising2D(L, mod.T, init=mode, seed=int(np.random.randint(1e6))).s
        etat["hist"].clear()

    sT.on_changed(lambda v: setattr(mod, "T", float(v)))
    sB.on_changed(lambda v: setattr(mod, "B", float(v)))
    radio.on_clicked(lambda lab: etat.__setitem__("algo", lab))
    b_alea.on_clicked(lambda _e: reinit("chaud"))
    b_ordre.on_clicked(lambda _e: reinit("froid"))
    b_neel.on_clicked(lambda _e: reinit("neel"))
    b_pause.on_clicked(lambda _e: etat.__setitem__("pause", not etat["pause"]))

    def bascule_J(_e=None):
        mod.J = -mod.J
        etat["hist"].clear()
        ax_m.set_ylabel("m" if mod.J > 0 else "m alterné")
    b_J.on_clicked(bascule_J)

    def maj(_k):
        if not etat["pause"]:
            algo = etat["algo"] if (mod.B == 0 and mod.J > 0) else "metropolis"
            mod.balayage(int(sV.val), algo)
        m = mod.aimantation() if mod.J > 0 else mod.aimantation_alternee()
        etat["hist"].append((mod.T, m))
        etat["hist"] = etat["hist"][-400:]
        im.set_data(mod.s)
        h = np.array(etat["hist"])
        pt.set_data(h[:, 0], h[:, 1]); cur.set_data([mod.T], [m])
        titre.set_text(f"T={mod.T:.2f}  B={mod.B:+.2f}  J={mod.J:+.0f}  {'m' if mod.J > 0 else 'm_alt'}={m:+.3f}  e={mod.energie():+.3f}  [{etat['algo']}]")
        return im, pt, cur, titre

    ani = FuncAnimation(fig, maj, interval=60, blit=False, cache_frame_data=False)
    fig._widgets = (sT, sB, sV, radio, b_alea, b_ordre, b_neel, b_pause, b_J)   # garde les références
    return fig, ani, mod


def lancer(**kw):
    fig, ani, _ = construire(**kw)
    plt.show()


# ----------------------------------------------------------------------------
# Démonstration scriptée -> GIF (sans écran)
# ----------------------------------------------------------------------------
SCENARIO = [
    # (légende, nb d'images, actions)  actions : dict(T=, B=, algo=, init=, vitesse=)
    ("1. Paramagnétique : T=3.5, spins désordonnés", 12, dict(T=3.5, B=0.0, init="chaud", vitesse=2)),
    ("2. Trempe à T=1.2 : les domaines de Weiss grossissent", 40, dict(T=1.2, vitesse=2)),
    ("3. On réchauffe vers Tc=2.269 : l'ordre fond", 30, dict(T=2.1, vitesse=3)),
    ("4. Au voisinage de Tc : amas à toutes les échelles (Wolff)", 22, dict(T=2.27, algo="wolff", vitesse=2)),
    ("5. T=3.5 : retour au désordre", 14, dict(T=3.5, algo="metropolis", vitesse=2)),
    ("6. Un champ B=+0.6 aligne les spins (T=1.8)", 26, dict(T=1.8, B=0.6, vitesse=2)),
    ("7. On inverse le champ : B=-0.6, nucléation de domaines ↓", 34, dict(B=-0.6, vitesse=2)),
    ("8. Antiferro (J→−J) : on part du damier de Néel, T=0.8, B=0", 14, dict(B=0.0, T=0.8, init="neel", J=-1, vitesse=1)),
    ("9. Antiferro chauffé à T=3 : l'ordre alterné disparaît", 22, dict(T=3.0, vitesse=2)),
]


def demo_gif(chemin="figures/15_demo_simulateur.gif", L=64, fps=10, scenario=SCENARIO, seed=1):
    """Pilote le simulateur interactif (curseurs, boutons, choix d'algo) et enregistre un GIF."""
    import io
    from pathlib import Path

    from PIL import Image

    fig, ani, mod = construire(L=L, seed=seed)
    sT, sB, sV, radio, b_alea, b_ordre, b_neel, _, b_J = fig._widgets
    fig.set_size_inches(10, 5.6)
    legende = fig.text(0.5, 0.975, "", ha="center", fontsize=13, weight="bold")
    boutons = {"chaud": b_alea, "froid": b_ordre, "neel": b_neel}
    images = []
    for texte, n, act in scenario:
        if "init" in act:
            boutons[act["init"]]._observers.process("clicked", None)
        if act.get("J") and act["J"] != mod.J:
            b_J._observers.process("clicked", None)
        if "T" in act:
            sT.set_val(act["T"])
        if "B" in act:
            sB.set_val(act["B"])
        if "vitesse" in act:
            sV.set_val(act["vitesse"])
        if "algo" in act:
            radio.set_active(0 if act["algo"] == "metropolis" else 1)
        legende.set_text(texte)
        for k in range(n):
            ani._func(k)
            fig.canvas.draw()
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=70)
            buf.seek(0)
            images.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=64))
    chemin = Path(chemin)
    chemin.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(chemin, save_all=True, append_images=images[1:], duration=int(1000 / fps), loop=0,
                   optimize=True)
    plt.close(fig)
    return chemin
