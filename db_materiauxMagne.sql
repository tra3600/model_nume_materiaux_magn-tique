/*6. Requête pour obtenir le nom de tous les matériaux avec une température de Curie strictement inférieure à 500 kelvins*/

SELECT nom 
FROM materiaux 
WHERE t_curie < 500;

/*7. Requête pour obtenir les noms de tous les fournisseurs proposant du nickel et le prix proposé par chacun pour 4,5 kilogrammes de nickel */

SELECT f.nom_fournisseur, p.prix_kg * 4.5 AS prix_total
FROM fournisseurs f
JOIN prix p ON f.id_fournisseur = p.id_four
WHERE p.id_mat = 8713;

/*8. Requête pour obtenir le nom du fournisseur de nickel le moins cher et le prix à payer pour 4,5 kilogrammes de nickel */

SELECT f.nom_fournisseur, p.prix_kg * 4.5 AS prix_total
FROM fournisseurs f
JOIN prix p ON f.id_fournisseur = p.id_four
WHERE p.id_mat = 8713
ORDER BY prix_total ASC
LIMIT 1;

/*9. Requête pour obtenir le nom de tous les matériaux et le prix moyen pour un kilogramme de chacun de ces matériaux avec une moyenne strictement inférieure à 50 euros par kilogramme
SQL*/

SELECT m.nom, AVG(p.prix_kg) AS prix_moyen
FROM materiaux m
JOIN prix p ON m.id_materiau = p.id_mat
GROUP BY m.nom
HAVING AVG(p.prix_kg) < 50;


