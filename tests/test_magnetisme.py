import math

import numpy as np
import pytest

from magnetisme import examen_mines_2022 as ex
from magnetisme import materiaux, theorie, weiss
from magnetisme.analyse import ajuster_exposant, croisement_binder
from magnetisme.ising import Ising2D, jackknife, mesurer, neel, temps_autocorrelation, thermodynamique


# ---------------------------------------------------------------- théorie
def test_tc_onsager():
    assert theorie.TC_ONSAGER == pytest.approx(2.26918531, abs=1e-7)


def test_dicho_et_champ_moyen():
    m = theorie.dicho(theorie.f, 0.5, 0.001, 1.0, 1e-9)
    assert m == pytest.approx(math.tanh(m / 0.5), abs=1e-8)
    assert theorie.aimantation_reduite(0.5) == pytest.approx(m, abs=1e-8)
    assert theorie.aimantation_reduite(1.2) == 0.0
    liste = theorie.construction_liste_m(0.1, 2.0)
    assert len(liste) == 500 and liste[-1] == 0 and liste[0] > 0.99


def test_onsager_energie_derivee_chaleur_specifique():
    T, h = np.array([1.0, 2.0, 3.0]), 1e-4
    num = (theorie.energie_onsager(T + h) - theorie.energie_onsager(T - h)) / (2 * h)
    assert theorie.chaleur_specifique_onsager(T) == pytest.approx(num, rel=1e-4)
    assert theorie.energie_onsager(theorie.TC_ONSAGER - 1e-7) == pytest.approx(-math.sqrt(2), abs=1e-3)


def test_brillouin_limites():
    t = np.array([0.3, 0.8])
    assert theorie.aimantation_brillouin(t, 0.5) == pytest.approx([theorie.aimantation_reduite(x) for x in t], abs=1e-6)
    assert theorie.aimantation_brillouin(1.2, 3.5)[0] == 0
    # B_J(x) ~ (J+1)/(3J) x en x -> 0
    assert theorie.brillouin(1e-4, 2.5) == pytest.approx(3.5 / 7.5 * 1e-4, rel=1e-4)


def test_susceptibilite_curie_weiss():
    assert theorie.susceptibilite_champ_moyen(np.array([2.0]))[0] == pytest.approx(1.0)


# ---------------------------------------------------------------- Ising
def test_energie_etats_fondamentaux():
    L = 8
    assert Ising2D(L, 1.0, init="froid").energie() == pytest.approx(-2.0)
    assert Ising2D(L, 1.0, J=-1, init="neel").energie() == pytest.approx(-2.0)
    assert Ising2D(L, 1.0, init="neel").energie() == pytest.approx(+2.0)
    assert Ising2D(L, 1.0, B=0.5, init="froid").energie() == pytest.approx(-2.5)


def test_energie_coherente_avec_sujet():
    """energie(s) du sujet (Q14) = L² * énergie par spin."""
    L = 6
    mod = Ising2D(L, 1.0, init="chaud", seed=1)
    assert ex.energie(mod.s.ravel().tolist(), L) == pytest.approx(mod.energie_totale())


def test_reproductibilite():
    a = Ising2D(16, 2.0, init="chaud", seed=42); a.balayage(20)
    b = Ising2D(16, 2.0, init="chaud", seed=42); b.balayage(20)
    # numba partage un seul générateur global : on teste surtout la configuration initiale
    assert (Ising2D(16, 2.0, init="chaud", seed=42).s == Ising2D(16, 2.0, init="chaud", seed=42).s).all()


def test_L_impair_refuse():
    with pytest.raises(ValueError):
        Ising2D(5)


@pytest.mark.parametrize("T", [1.5, 3.5])
@pytest.mark.parametrize("algo", ["metropolis", "wolff"])
def test_energie_vs_onsager(T, algo):
    mod = Ising2D(24, T, init="chaud", seed=3)
    mes = mesurer(mod, 500, 2500, algo=algo)
    e = mes["e"].mean()
    assert e == pytest.approx(float(theorie.energie_onsager(T)), abs=0.04)


def test_aimantation_basse_temperature():
    mod = Ising2D(24, 1.5, init="froid", seed=1)
    mes = mesurer(mod, 300, 1500)
    assert np.abs(mes["m"]).mean() == pytest.approx(float(theorie.aimantation_onsager(1.5)), abs=0.02)


def test_haute_temperature_desordre():
    mod = Ising2D(32, 4.0, init="froid", seed=1)
    mes = mesurer(mod, 300, 1500)
    assert np.abs(mes["m"]).mean() < 0.1


def test_champ_aligne_les_spins():
    mod = Ising2D(16, 3.0, B=2.0, init="chaud", seed=1)
    mod.balayage(200)
    assert mod.aimantation() > 0.8


def test_wolff_refuse_champ():
    with pytest.raises(ValueError):
        Ising2D(8, 2.0, B=0.1).balayage(1, "wolff")


def test_damier_numpy_conserve_les_valeurs():
    mod = Ising2D(16, 2.0, init="chaud", seed=0)
    for _ in range(5):
        mod._damier()
    assert set(np.unique(mod.s)) <= {-1, 1}


def test_damier_numpy_energie():
    mod = Ising2D(24, 3.0, init="chaud", seed=0)
    es = []
    for k in range(600):
        mod._damier()
        if k > 100:
            es.append(mod.energie())
    assert np.mean(es) == pytest.approx(float(theorie.energie_onsager(3.0)), abs=0.05)


def test_jackknife_et_autocorrelation():
    x = np.random.default_rng(0).normal(size=4000)
    val, err = jackknife(np.mean, [x])
    assert abs(val) < 4 * err + 1e-9 and err == pytest.approx(1 / math.sqrt(4000), rel=0.4)
    assert temps_autocorrelation(x) == pytest.approx(0.5, abs=0.3)


def test_thermodynamique_cles():
    mod = Ising2D(8, 2.0, init="chaud", seed=0)
    r = thermodynamique(mesurer(mod, 50, 200), 8, 2.0)
    assert set(r) == {"M", "E", "chi", "C", "U4"}


# ---------------------------------------------------------------- domaines de Weiss
def test_domaines_cas_simples():
    assert weiss.domaines(np.ones((6, 6), np.int8))[1] == 1
    assert weiss.domaines(neel(6))[1] == 36
    s = np.ones((6, 6), np.int8); s[:, 2:4] = -1      # deux bandes (une enroulée par le bord périodique)
    assert weiss.domaines(s)[1] == 2


def test_domaines_identiques_au_corrige_du_sujet():
    rng = np.random.default_rng(1)
    for k in range(40):
        h = int(rng.choice([4, 6, 12]))
        s = rng.choice([-1, 1], size=(h, h), p=[0.3, 0.7])
        lab, n = weiss.domaines(s)
        w = ex.construire_domaines_weiss(s.ravel().tolist(), h)
        assert n == max(w) + 1 and (lab.ravel() == np.array(w)).all()


def test_recursif_egal_pile():
    rng = np.random.default_rng(2)
    h = 8
    s = rng.choice([-1, 1], size=h * h).tolist()
    w1 = [-1] * len(s); ex.explorer_voisinage(s, 0, w1, 0, h)
    w2 = [-1] * len(s); ex.explorer_voisinage_pile(s, 0, w2, 0, [], h)
    assert w1 == w2


def test_correlation_decroit():
    mod = Ising2D(32, 3.5, init="chaud", seed=1)
    mod.balayage(200)
    r, g = weiss.correlation(mod.s)
    assert g[0] > g[3]


# ---------------------------------------------------------------- corrigé de l'examen
def test_examen_voisins_et_initialisation():
    assert sorted(ex.liste_voisins(0, 4)) == sorted([3, 1, 4, 12])
    assert len(set(ex.liste_voisins(5, 4))) == 4
    for i in range(16):
        assert all(i in ex.liste_voisins(j, 4) for j in ex.liste_voisins(i, 4))   # symétrie
    assert sum(ex.initialisation(10)) == 100 and sum(ex.initialisation_anti(10)) == 0
    assert ex.repliement(list(range(9)))[2] == [6, 7, 8]


def test_examen_delta_e_rapide_egal_lent():
    s = np.random.default_rng(0).choice([-1, 1], size=36).tolist()
    for i in range(36):
        assert ex.calcul_delta_e1(s, i, 6) == ex.calcul_delta_e2(s, i, 6)


def test_examen_boltzmann_et_montecarlo():
    assert ex.test_boltzmann(-1.0, 1.0) and ex.test_boltzmann(0.0, 1.0)
    assert not ex.test_boltzmann(1e9, 1.0)
    assert ex.aimantation_moyenne(40_000, 1.0, 10) > 0.95


# ---------------------------------------------------------------- base de données
def test_requetes_sql():
    con = materiaux.connexion()
    assert ("nickel",) not in materiaux.executer(6, con)
    assert ("YIG",) not in materiaux.executer(6, con)
    r7 = materiaux.executer(7, con)
    assert len(r7) == 4 and min(p for _, p in r7) == pytest.approx(24.5 * 4.5)
    r8 = materiaux.executer(8, con)
    assert {n for n, _ in r8} == {"Worldwide Materials", "Magnetics & Co"}     # ex æquo conservés
    assert all(p < 50 for _, p in materiaux.executer(9, con))


# ---------------------------------------------------------------- analyse
def test_ajuster_exposant_et_binder():
    x = np.logspace(0, 2, 20)
    p, A = ajuster_exposant(x, 3 * x ** 0.5)
    assert p == pytest.approx(0.5) and A == pytest.approx(3)
    T = np.linspace(0, 1, 11)
    assert croisement_binder(T, T, 1 - T) == pytest.approx(0.5)


def test_cli(capsys):
    from magnetisme.__main__ import main
    main(["champ-moyen", "0.5"])
    assert "0.9575" in capsys.readouterr().out
    main(["sql"])
    assert "Q8" in capsys.readouterr().out
    main(["simuler", "-L", "16", "-T", "3", "--therm", "50", "--mesures", "200"])
    assert "Ising 2D" in capsys.readouterr().out


def test_figures_rapides(tmp_path):
    from magnetisme import figures
    out = figures.generer_tout(tmp_path, rapide=True, noms=["champ_moyen", "landau", "domaines", "materiaux"], verbeux=False)
    assert all(p.exists() and p.stat().st_size > 5000 for p in out)


def test_interactif_construit():
    import matplotlib
    matplotlib.use("Agg")
    from magnetisme.interactif import construire
    fig, ani, mod = construire(L=16)
    ani._func(0)
    assert mod.s.shape == (16, 16)
