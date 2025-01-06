from math import exp, tanh
from random import randrange, random

def f(x, t):
    """
    Calcule la valeur de l'équation f(x, t) = tanh(x/t) - x.
    
    Paramètres:
    x (float): Variable inconnue.
    t (float): Température réduite.
    
    Retourne:
    float: Valeur de f(x, t).
    """
    return tanh(x / t) - x

def dicho(f, t, a, b, eps):
    """
    Calcule une valeur approchée du zéro d'une fonction f(m, t) sur l'intervalle [a, b] à eps près.
    
    Paramètres:
    f (callable): Fonction dont on cherche le zéro.
    t (float): Paramètre de la fonction f.
    a (float): Borne inférieure de l'intervalle de recherche.
    b (float): Borne supérieure de l'intervalle de recherche.
    eps (float): Précision de la recherche.
    
    Retourne:
    float: Valeur approchée du zéro de la fonction f.
    """
    while (b - a) / 2.0 > eps:
        midpoint = (a + b) / 2.0
        if f(midpoint, t) == 0:
            return midpoint
        elif f(a, t) * f(midpoint, t) < 0:
            b = midpoint
        else:
            a = midpoint
    return (a + b) / 2.0

# Exemple d'utilisation
t = 0.5
a = 0.001
b = 1.0
eps = 1e-6
m_approx = dicho(f, t, a, b, eps)
print(m_approx)


def construction_liste_m(t1, t2):
    """
    Construit une liste de 500 solutions de l'équation (1) pour t variant linéairement de t1 à t2.
    
    Paramètres:
    t1 (float): Valeur initiale de t.
    t2 (float): Valeur finale de t.
    
    Retourne:
    List[float]: Liste des solutions m pour chaque valeur de t.
    """
    m_list = []
    t_values = [t1 + i * (t2 - t1) / 499 for i in range(500)]
    for t in t_values:
        if t > 1:
            m_list.append(0)
        else:
            m_approx = dicho(f, t, 0.001, 1.0, 1e-6)
            m_list.append(m_approx)
    return m_list

# Exemple d'utilisation
t1 = 0.1
t2 = 2.0
m_list = construction_liste_m(t1, t2)
print(m_list)

