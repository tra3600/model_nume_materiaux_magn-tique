"""Interface en ligne de commande :  python -m magnetisme <commande> [options]"""
from __future__ import annotations

import argparse
import sys


def main(argv=None):
    p = argparse.ArgumentParser(prog="magnetisme", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("figures", help="génère toutes les figures dans figures/")
    s.add_argument("--rapide", action="store_true", help="tailles réduites (~15 s)")
    s.add_argument("--dossier", default="figures")
    s.add_argument("--only", nargs="*", help="sous-ensemble : " + ", ".join(
        __import__("magnetisme.figures", fromlist=["TOUTES"]).TOUTES))

    s = sub.add_parser("simuler", help="une simulation Ising et ses observables")
    s.add_argument("-L", type=int, default=64)
    s.add_argument("-T", type=float, default=2.0)
    s.add_argument("-B", type=float, default=0.0)
    s.add_argument("--J", type=float, default=1.0)
    s.add_argument("--algo", choices=["metropolis", "wolff"], default="metropolis")
    s.add_argument("--therm", type=int, default=1000)
    s.add_argument("--mesures", type=int, default=5000)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--image", help="enregistre la configuration finale (PNG)")

    s = sub.add_parser("balayage", help="tableau M, E, chi, C, U4 vs T (CSV sur stdout)")
    s.add_argument("-L", type=int, default=32)
    s.add_argument("--tmin", type=float, default=1.5)
    s.add_argument("--tmax", type=float, default=3.5)
    s.add_argument("--n", type=int, default=21)
    s.add_argument("--mesures", type=int, default=5000)

    s = sub.add_parser("domaines", help="domaines de Weiss d'une configuration")
    s.add_argument("-L", type=int, default=100)
    s.add_argument("-T", type=float, default=2.269)
    s.add_argument("--image")

    sub.add_parser("sql", help="requêtes 6 à 9 du sujet sur la base de matériaux")
    s = sub.add_parser("champ-moyen", help="résout m = tanh(m/t) (partie I du sujet)")
    s.add_argument("t", type=float, nargs="+")
    sub.add_parser("examen", help="exécute le corrigé Python pur du sujet Mines 2022")
    sub.add_parser("interactif", help="simulateur interactif (curseurs T, B)")
    s = sub.add_parser("demo", help="GIF de démonstration du simulateur (sans écran)")
    s.add_argument("--sortie", default="figures/15_demo_simulateur.gif")
    s.add_argument("-L", type=int, default=64)
    a = p.parse_args(argv)

    if a.cmd == "figures":
        from .figures import generer_tout
        generer_tout(a.dossier, rapide=a.rapide, noms=a.only)
    elif a.cmd == "simuler":
        import numpy as np
        from .ising import HAS_NUMBA, Ising2D, mesurer, thermodynamique
        from .theorie import aimantation_onsager, energie_onsager
        mod = Ising2D(a.L, a.T, J=a.J, B=a.B, init="chaud", seed=a.seed)
        mes = mesurer(mod, a.therm, a.mesures, algo=a.algo)
        r = thermodynamique(mes, a.L, a.T)
        print(f"Ising 2D  L={a.L}  T={a.T}  J={a.J}  B={a.B}  [{a.algo}, numba={'oui' if HAS_NUMBA else 'non'}]")
        for k, nom in (("M", "<|m|>"), ("E", "énergie/spin"), ("chi", "susceptibilité"),
                       ("C", "chaleur spéc."), ("U4", "Binder U4")):
            print(f"  {nom:15s} = {r[k][0]: .5f} ± {r[k][1]:.5f}")
        if a.B == 0 and a.J == 1:
            print(f"  exact (L=∞)     : |m| = {float(aimantation_onsager(a.T)):.5f}, e = {float(energie_onsager(a.T)):.5f}")
        if a.image:
            import matplotlib; matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from .figures import SPINS
            plt.imsave(a.image, mod.s, cmap=SPINS, vmin=-1, vmax=1)
            print("image :", a.image)
    elif a.cmd == "balayage":
        import numpy as np
        from .analyse import balayage_temperature
        T = np.linspace(a.tmin, a.tmax, a.n)
        r = balayage_temperature(a.L, T, a.mesures // 5, a.mesures)
        print("T,M,errM,E,errE,chi,errchi,C,errC,U4,errU4")
        for k, t in enumerate(T):
            print(f"{t:.4f}," + ",".join(f"{r[g][k, 0]:.5f},{r[g][k, 1]:.5f}" for g in ("M", "E", "chi", "C", "U4")))
    elif a.cmd == "domaines":
        from .ising import Ising2D
        from .weiss import statistiques, domaines
        mod = Ising2D(a.L, a.T, init="chaud", seed=0)
        mod.balayage(2000, "wolff")
        for k, v in statistiques(mod.s).items():
            print(f"  {k:16s} {v:.4g}" if isinstance(v, float) else f"  {k:16s} {v}")
        if a.image:
            import matplotlib; matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
            lab, n = domaines(mod.s)
            plt.imsave(a.image, np.random.default_rng(0).permutation(n)[lab], cmap="nipy_spectral")
            print("image :", a.image)
    elif a.cmd == "sql":
        from .materiaux import REQUETES, connexion
        con = connexion()
        for q, (titre, sql) in REQUETES.items():
            print(f"\nQ{q}. {titre}")
            for row in con.execute(sql):
                print("   ", row)
    elif a.cmd == "champ-moyen":
        from .theorie import aimantation_reduite
        for t in a.t:
            print(f"t = {t:g}  ->  m = {aimantation_reduite(t):.6f}")
    elif a.cmd == "examen":
        from . import examen_mines_2022
        import runpy
        runpy.run_module("magnetisme.examen_mines_2022", run_name="__main__")
    elif a.cmd == "demo":
        import matplotlib; matplotlib.use("Agg")
        from .interactif import demo_gif
        print("GIF :", demo_gif(a.sortie, L=a.L))
    elif a.cmd == "interactif":
        from .interactif import lancer
        lancer()
    return 0


if __name__ == "__main__":
    sys.exit(main())
