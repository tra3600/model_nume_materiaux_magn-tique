-- Températures de Curie / Néel : valeurs de la littérature (arrondies). Prix : FICTIFS.
INSERT INTO materiaux VALUES
 (4534, 'cobalt',            1388, 'ferromagnétique',  0.5, 1440),
 (1254, 'dioxyde de chrome',  386, 'ferromagnétique',  1.0,  480),
 (8713, 'nickel',             627, 'ferromagnétique',  0.5,  485),
 (8284, 'YIG',                560, 'ferrimagnétique',  2.5,  194),
 (1001, 'fer',               1043, 'ferromagnétique',  0.5, 1750),
 (1002, 'gadolinium',         292, 'ferromagnétique',  3.5, 2060),
 (1003, 'magnétite',          858, 'ferrimagnétique',  2.0,  480),
 (1004, 'oxyde d''europium',   69, 'ferromagnétique',  3.5, 1900),
 (1005, 'dysprosium',          88, 'ferromagnétique',  7.5, 2900),
 (1006, 'MnBi',               630, 'ferromagnétique',  1.0,  620),
 (1007, 'chrome',             311, 'antiferromagnétique', 1.5, 0),
 (1008, 'NiO',                523, 'antiferromagnétique', 1.0, 0);
INSERT INTO fournisseurs VALUES
 (145, 'Worldwide Materials'), (13, 'Materials Company'), (77, 'Magnetics & Co'), (42, 'EuroAlloy');
INSERT INTO prix VALUES
 (1, 8713, 145, 24.50), (2, 8713,  13, 27.90), (3, 8713, 77, 24.50), (4, 8713, 42, 31.00),
 (5, 4534, 145, 38.00), (6, 4534,  13, 41.20), (7, 1001, 145,  0.90), (8, 1001, 42,  1.10),
 (9, 1254,  77, 52.75), (10, 1254, 13, 50.40), (11, 8284, 42, 1357.30), (12, 1002, 145, 120.00),
 (13, 1003, 77,  1.80), (14, 1003, 13,  2.10), (15, 1006, 42, 95.00), (16, 1004, 77, 900.00),
 (17, 1005, 145, 350.00), (18, 1005, 13, 340.00);
