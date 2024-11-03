# Séances 4-5 : Adaptation de Maillage Instationnaire et Méthodes d’Intégration

Dans ces séances, l’objectif est d’étendre le modèle d’adaptation de maillage pour des solutions instationnaires en modifiant le code et en analysant l'impact des différents schémas temporels et du maillage sur l'erreur.

## **Code `adrs_insta.py` : Solution Instationnaire avec Solution Exacte Cible Instationnaire**
- Visualiser l’erreur au milieu et à la fin du calcul (temps T/2 et T) pour différents maillages uniformes.
- Analyser l’évolution de l’erreur au point milieu du domaine pour divers schémas de Runge-Kutta (ordre 1 à 4).

## **Code `adrs_insta_multiple_mesh_adap.py` : Adaptation de Maillage avec Terme Source Dépendant du Temps**
- Modification du code pour inclure un terme source temporel, défini par une solution exacte \( u(t, s) = u(t) v(s) \), avec :
  - \( u(t)v(s) = \sin(4\pi t) v(s) \)
  - \( u'v = 4\pi \cos(4\pi t) v(s) \)
- Visualisation de la solution à divers instants pour une période de 1 seconde.
- Introduction de critères d'arrêt mixtes basés sur le nombre de points de maillage et l'erreur \( L2 \).
- Analyse du maillage adaptatif en utilisant une métrique basée sur la solution finale à \( t = \text{Time} \).
- Inclusion de l’adaptation instationnaire avec une moyenne temporelle des métriques pour optimiser le maillage.

---

