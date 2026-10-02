"""Base SQLite de matériaux magnétiques et requêtes de la partie II du sujet (questions 6 à 9)."""
from __future__ import annotations

import sqlite3
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"

REQUETES = {
    6: ("Matériaux de température de Curie < 500 K",
        "SELECT nom FROM materiaux WHERE t_curie < 500 ORDER BY nom"),
    7: ("Fournisseurs de nickel et prix de 4,5 kg",
        """SELECT f.nom_fournisseur, ROUND(p.prix_kg * 4.5, 2) AS prix_total
           FROM fournisseurs f JOIN prix p ON f.id_fournisseur = p.id_four
           WHERE p.id_mat = 8713 ORDER BY prix_total"""),
    8: ("Fournisseur(s) de nickel le(s) moins cher(s) (ex æquo conservés)",
        """SELECT f.nom_fournisseur, ROUND(p.prix_kg * 4.5, 2) AS prix_total
           FROM fournisseurs f JOIN prix p ON f.id_fournisseur = p.id_four
           WHERE p.id_mat = 8713
             AND p.prix_kg = (SELECT MIN(prix_kg) FROM prix WHERE id_mat = 8713)"""),
    9: ("Prix moyen au kg < 50 € par matériau",
        """SELECT m.nom, ROUND(AVG(p.prix_kg), 2) AS prix_moyen
           FROM materiaux m JOIN prix p ON m.id_materiau = p.id_mat
           GROUP BY m.id_materiau HAVING AVG(p.prix_kg) < 50 ORDER BY prix_moyen"""),
}


def connexion() -> sqlite3.Connection:
    """Base en mémoire construite à partir de data/schema.sql et data/donnees.sql."""
    con = sqlite3.connect(":memory:")
    con.executescript((DATA / "schema.sql").read_text(encoding="utf-8"))
    con.executescript((DATA / "donnees.sql").read_text(encoding="utf-8"))
    return con


def executer(numero: int, con: sqlite3.Connection | None = None):
    con = con or connexion()
    return con.execute(REQUETES[numero][1]).fetchall()


def table_materiaux(con: sqlite3.Connection | None = None):
    """Liste de dicts (nom, t_curie, ordre, j_eff, ms_kAm, prix_moyen) pour les tracés."""
    con = con or connexion()
    cur = con.execute("""SELECT m.nom, m.t_curie, m.ordre, m.j_eff, m.ms_kAm, AVG(p.prix_kg)
                         FROM materiaux m LEFT JOIN prix p ON p.id_mat = m.id_materiau
                         GROUP BY m.id_materiau ORDER BY m.t_curie DESC""")
    cles = ("nom", "t_curie", "ordre", "j_eff", "ms_kAm", "prix_moyen")
    return [dict(zip(cles, r)) for r in cur]
