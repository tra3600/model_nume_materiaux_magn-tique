def initialisation() -> list:
    """
    Renvoie une liste d'initialisation des domaines contenant n spins de valeur 1.
    
    Retourne:
    list: Liste de n éléments, chacun valant 1.
    """
    h = 100
    n = h ** 2
    return [1] * n

# Exemple d'utilisation
spins = initialisation()
print(spins)


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

# Exemple d'utilisation
anti_spins = initialisation_anti()
print(anti_spins)


def repliement(s: list) -> list:
    """
    Convertit une liste de spins en un tableau de taille h x h.
    
    Paramètres:
    s (list): Liste de spins.
    
    Retourne:
    list: Liste de h listes de taille h représentant le domaine.
    """
    h = int(len(s) ** 0.5)
    return [s[i * h:(i + 1) * h] for i in range(h)]

# Exemple d'utilisation
replie_spins = repliement(initialisation())
print(replie_spins)


def liste_voisins(i: int, h: int) -> list:
    """
    Renvoie la liste des indices des plus proches voisins du spin s_i d'indice i dans la liste s.
    
    Paramètres:
    i (int): Indice du spin.
    h (int): Taille de la grille (nombre de spins par ligne/colonne).
    
    Retourne:
    list: Liste des indices des voisins (gauche, droite, dessous, dessus).
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

# Exemple d'utilisation
i = 5
h = 10
voisins = liste_voisins(i, h)
print(voisins)


def energie(s: list, h: int) -> float:
    """
    Calcule l'énergie d'une configuration s donnée.
    
    Paramètres:
    s (list): Liste de spins.
    h (int): Taille de la grille (nombre de spins par ligne/colonne).
    
    Retourne:
    float: Énergie de la configuration.
    """
    E = 0
    for i in range(len(s)):
        for j in liste_voisins(i, h):
            E -= s[i] * s[j]
    return E / 2  # Diviser par 2 pour éviter de compter deux fois chaque interaction

# Exemple d'utilisation
h = 10
spins = initialisation()
E = energie(spins, h)
print(E)

import random
import math

def test_boltzmann(delta_e: float, T: float) -> bool:
    """
    Teste si un spin change de signe en utilisant la loi de Boltzmann.
    
    Paramètres:
    delta_e (float): Variation d'énergie.
    T (float): Température.
    
    Retourne:
    bool: True si le spin change de signe, False sinon.
    """
    if delta_e <= 0:
        return True
    else:
        p = math.exp(-delta_e / T)
        return random.random() < p

# Exemple d'utilisation
delta_e = 0.5
T = 300
print(test_boltzmann(delta_e, T))


def calcul_delta_e2(s: list, i: int, h: int) -> float:
    """
    Calcule la variation d'énergie delta_e résultant d'un changement d'orientation du spin i.
    
    Paramètres:
    s (list): Liste des spins.
    i (int): Indice du spin à changer.
    h (int): Taille de la grille (nombre de spins par ligne/colonne).
    
    Retourne:
    float: Variation d'énergie delta_e.
    """
    delta_e = 0
    for j in liste_voisins(i, h):
        delta_e += 2 * s[i] * s[j]
    return delta_e

# Exemple d'utilisation
h = 10
spins = initialisation()
i = random.randint(0, len(spins) - 1)
delta_e = calcul_delta_e2(spins, i, h)
print(delta_e)


def monte_carlo(s: list, T: float, n_tests: int, h: int) -> list:
    """
    Applique la méthode de Monte-Carlo pour faire évoluer la liste des spins s.
    
    Paramètres:
    s (list): Liste des spins.
    T (float): Température.
    n_tests (int): Nombre de tests de Boltzmann à effectuer.
    h (int): Taille de la grille (nombre de spins par ligne/colonne).
    
    Retourne:
    list: Liste des spins modifiée.
    """
    for _ in range(n_tests):
        i = random.randint(0, len(s) - 1)
        delta_e = calcul_delta_e2(s, i, h)
        if test_boltzmann(delta_e, T):
            s[i] = -s[i]
    return s

# Exemple d'utilisation
T = 300
n_tests = 1000
spins = initialisation()
h = 10
spins = monte_carlo(spins, T, n_tests, h)
print(spins)



def aimantation_moyenne(n_tests: int, T: float) -> float:
    """
    Calcule l'aimantation moyenne d'une configuration de spins à la température T.
    
    Paramètres:
    n_tests (int): Nombre de tests de Boltzmann à effectuer.
    T (float): Température.
    
    Retourne:
    float: Aimantation moyenne de la configuration.
    """
    spins = initialisation()
    h = 100
    spins = monte_carlo(spins, T, n_tests, h)
    aimantation = sum(spins) / len(spins)
    return aimantation

# Exemple d'utilisation
T = 300
n_tests = 1000
aimantation = aimantation_moyenne(n_tests, T)
print(aimantation)



