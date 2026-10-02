"""Corrigé commenté de l'épreuve d'informatique Mines-Ponts 2022 (MP-PC-PSI) en Python pur.

Sujet : « Modélisation numérique d'un matériau magnétique » (data/PSI_INFO_MINES_1_2022.enonce.pdf).
Contraintes du sujet respectées : pas de numpy, spins stockés dans une liste 1D de n = h² entiers.
Aucun effet de bord à l'import (les anciennes versions affichaient des listes de 10 000 éléments).

Chaque fonction porte le numéro de la question correspondante.
"""
from math import exp, tanh            # Q1
from random import randrange, random  # Q1

h = 100          # côté de l'échantillon (variable globale du sujet)
n = h ** 2


# --------------------------------------------------------------------------- Partie I
def f(x, t):                         # Q2 : m = tanh(m/t)  <=>  f(m,t) = tanh(m/t) - m = 0
    return tanh(x / t) - x


def dicho(f, t, a, b, eps):          # Q3
    """Dichotomie. Complexité (Q4) : Θ(log2((b-a)/eps)) évaluations de f."""
    fa = f(a, t)
    while (b - a) / 2 > eps:
        m = (a + b) / 2
        fm = f(m, t)
        if fm == 0:
            return m
        if fa * fm < 0:
            b = m
        else:
            a, fa = m, fm
    return (a + b) / 2


def construction_liste_m(t1, t2):    # Q5
    ts = [t1 + k * (t2 - t1) / 499 for k in range(500)]   # 500 valeurs, bornes incluses
    return [0 if t >= 1 else dicho(f, t, 0.001, 1, 1e-6) for t in ts]


# --------------------------------------------------------------------------- Partie III
def initialisation(h=h):             # Q10 : tous les spins up
    return [1] * h ** 2


def initialisation_anti(h=h):        # Q11 : damier +1/-1 (h pair)
    return [1 if (i // h + i % h) % 2 == 0 else -1 for i in range(h * h)]


def repliement(s):                   # Q12 : liste 1D -> liste de h listes (affichage uniquement)
    h = int(len(s) ** 0.5)
    return [s[i * h:(i + 1) * h] for i in range(h)]


def liste_voisins(i, h=h):           # Q13 : gauche, droite, dessous, dessus (bords périodiques)
    l, c = i // h, i % h
    return [l * h + (c - 1) % h, l * h + (c + 1) % h, ((l + 1) % h) * h + c, ((l - 1) % h) * h + c]


def energie(s, h=h):                 # Q14 : E = -(J/2) sum_i sum_{j in V_i} s_i s_j, J = 1
    return -sum(s[i] * s[j] for i in range(len(s)) for j in liste_voisins(i, h)) / 2


def test_boltzmann(delta_e, T):      # Q15
    return delta_e <= 0 or random() < exp(-delta_e / T)


def calcul_delta_e1(s, i, h=h):      # Q16 (solution lente) : O(n) car recalcule toute l'énergie
    s2 = s[:]
    s2[i] = -s[i]
    return energie(s2, h) - energie(s, h)


def calcul_delta_e2(s, i, h=h):      # Q16 (solution retenue) : O(1), seuls 4 voisins comptent
    return sum(2 * s[i] * s[j] for j in liste_voisins(i, h))


def monte_carlo(s, T, n_tests, h=h):  # Q17 : modifie s en place
    for _ in range(n_tests):
        i = randrange(len(s))
        if test_boltzmann(calcul_delta_e2(s, i, h), T):
            s[i] = -s[i]


def aimantation_moyenne(n_tests, T, h=h):   # Q18 ; complexité Q19 : Θ(n + n_tests)
    s = initialisation(h)
    monte_carlo(s, T, n_tests, h)
    return sum(s) / len(s)


# Q20 : avec toutes les paires, delta_e coûterait Θ(n) au lieu de Θ(1) -> Θ(n_tests * n).
# Q21 : quand T croît, l'agitation thermique désordonne les spins : l'aimantation décroît, les
#       domaines se fragmentent, et au-delà de Tc ≈ 2.269 l'aimantation moyenne s'annule
#       (paramagnétisme).


# --------------------------------------------------------------------------- Partie IV
def explorer_voisinage(s, i, weiss, num, h=h):      # Q22 : version récursive
    weiss[i] = num
    for j in liste_voisins(i, h):
        if s[j] == s[i] and weiss[j] == -1:
            explorer_voisinage(s, j, weiss, num, h)     # attention : profondeur de récursion ~ taille du domaine


def explorer_voisinage_pile(s, i, weiss, num, pile, h=h):   # Q23 : pile explicite
    pile.append(i)
    weiss[i] = num
    while pile:
        k = pile.pop()
        for j in liste_voisins(k, h):
            if s[j] == s[k] and weiss[j] == -1:
                weiss[j] = num              # on marque à l'empilement : jamais empilé deux fois
                pile.append(j)


def construire_domaines_weiss(s, h=None):    # Q24
    h = h or int(len(s) ** 0.5)
    weiss = [-1] * len(s)
    num = 0
    for i in range(len(s)):
        if weiss[i] == -1:
            explorer_voisinage_pile(s, i, weiss, num, [], h)
            num += 1
    return weiss


if __name__ == "__main__":
    print("m(t=0.5)       =", dicho(f, 0.5, 0.001, 1, 1e-6))
    s = initialisation_anti(10)
    print("domaines (anti) =", max(construire_domaines_weiss(s, 10)) + 1)
    print("m(T=1.5, 20x20) =", aimantation_moyenne(200_000, 1.5, 20))
