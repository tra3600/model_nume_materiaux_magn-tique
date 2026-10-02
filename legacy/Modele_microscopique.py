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

