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


def initialisation_anti() -> list:
    """
    Renvoie une liste d'initialisation des domaines contenant h spins en largeur et h en hauteur en alternant les 1 et -1.
    
    Retourne:
    list: Liste de n éléments, alternant les 1 et -1.
    """
    h = 100
    n = h ** 2
    s = []
    for i in range(h):
        for j in range(h):
            if (i + j) % 2 == 0:
                s.append(1)
            else:
                s.append(-1)
    return s



def liste_voisins(i: int, h: int) -> list:
    """
    Renvoie la liste des indices des plus proches voisins du spin s_i d'indice i dans la liste s.
    """
    voisins = []
    # gauche
    voisins.append(i - 1 if i % h != 0 else i + h - 1)
    # droite
    voisins.append(i + 1 if (i + 1) % h != 0 else i - h + 1)
    # dessous
    voisins.append(i + h if i + h < h * h else i % h)
    # dessus
    voisins.append(i - h if i - h >= 0 else (h * (h - 1)) + (i % h))
    return voisins

def explorer_voisinage(s: list, i: int, weiss: list, num: int, h: int):
    """
    Explore récursivement le voisinage d'un spin pour construire la liste weiss.
    
    Paramètres:
    s (list): Liste des spins.
    i (int): Indice du spin de départ.
    weiss (list): Liste des domaines de Weiss.
    num (int): Numéro du domaine de Weiss.
    h (int): Taille de la grille (nombre de spins par ligne/colonne).
    """
    # Marquer le spin courant avec le numéro du domaine
    weiss[i] = num
    
    # Explorer les voisins
    for voisin in liste_voisins(i, h):
        # Vérifier si le voisin a la même valeur et n'a pas encore été affecté à un domaine
        if s[voisin] == s[i] and weiss[voisin] == -1:
            explorer_voisinage(s, voisin, weiss, num, h)

# Exemple d'utilisation
h = 10
s = initialisation()
weiss = [-1] * (h * h)
explorer_voisinage(s, 0, weiss, 0, h)
print(weiss)


def explorer_voisinage_pile(s: list, i: int, weiss: list, num: int, h: int, pile: list):
    """
    Explore itérativement le voisinage d'un spin pour construire la liste weiss en utilisant une pile explicite.
    
    Paramètres:
    s (list): Liste des spins.
    i (int): Indice du spin de départ.
    weiss (list): Liste des domaines de Weiss.
    num (int): Numéro du domaine de Weiss.
    h (int): Taille de la grille (nombre de spins par ligne/colonne).
    pile (list): Pile explicite pour l'exploration.
    """
    # Initialiser la pile avec le spin de départ
    pile.append(i)
    
    while pile:
        # Récupérer l'indice du spin à explorer
        current = pile.pop()
        # Marquer le spin courant avec le numéro du domaine
        weiss[current] = num
        
        # Explorer les voisins
        for voisin in liste_voisins(current, h):
            # Vérifier si le voisin a la même valeur et n'a pas encore été affecté à un domaine
            if s[voisin] == s[current] and weiss[voisin] == -1:
                pile.append(voisin)

# Exemple d'utilisation
h = 10
s = initialisation()
weiss = [-1] * (h * h)
pile = []
explorer_voisinage_pile(s, 0, weiss, 0, h, pile)
print(weiss)


def construire_domaines_weiss(s: list) -> list:
    """
    Construit et renvoie la liste weiss contenant le numéro des domaines de Weiss de chaque spin du domaine.
    
    Paramètres:
    s (list): Liste des spins.
    
    Retourne:
    list: Liste weiss contenant le numéro des domaines de Weiss de chaque spin.
    """
    h = int(len(s) ** 0.5)
    weiss = [-1] * len(s)
    num = 0
    
    for i in range(len(s)):
        if weiss[i] == -1:
            pile = []
            explorer_voisinage_pile(s, i, weiss, num, h, pile)
            num += 1
    
    return weiss

# Exemple d'utilisation
s = initialisation()
weiss = construire_domaines_weiss(s)
print(weiss)