-- Schéma de la base de matériaux magnétiques (partie II du sujet Mines 2022)
CREATE TABLE materiaux (
    id_materiau INTEGER PRIMARY KEY,
    nom         TEXT    NOT NULL,
    t_curie     INTEGER NOT NULL,          -- kelvin
    ordre       TEXT    NOT NULL DEFAULT 'ferromagnétique',
    j_eff       REAL,                      -- moment cinétique effectif (modèle de Brillouin)
    ms_kAm      REAL                       -- aimantation à saturation (kA/m, ≈ 0 K)
);
CREATE TABLE fournisseurs (
    id_fournisseur  INTEGER PRIMARY KEY,
    nom_fournisseur TEXT NOT NULL
);
CREATE TABLE prix (
    id_prix  INTEGER PRIMARY KEY,
    id_mat   INTEGER NOT NULL REFERENCES materiaux(id_materiau),
    id_four  INTEGER NOT NULL REFERENCES fournisseurs(id_fournisseur),
    prix_kg  REAL    NOT NULL              -- euros / kg (valeurs fictives)
);
