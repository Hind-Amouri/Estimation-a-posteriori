# Estimation a Posteriori et Adaptation de Maillage pour les EDP

Ce dépôt contient les codes Python développés dans le cadre du projet portant sur la modélisation numérique des équations aux dérivées partielles (EDP) avec estimation a posteriori de l'erreur et optimisation de l'adaptation de maillage.

# Compte Rendu des Séances de Cours

## Introduction
L'analyse numérique des équations aux dérivées partielles (EDP) est un domaine fondamental dans le champ des mathématiques appliquées. Ce cours explore les méthodes numériques permettant de résoudre ces équations complexes qui modélisent divers phénomènes physiques, tels que la diffusion, l'écoulement de fluides et la chaleur. L'objectif principal est de développer les compétences nécessaires pour aborder ces problèmes avec rigueur et précision.

## Motivation
Les EDP sont omniprésentes dans des disciplines telles que la physique, l'ingénierie, la biologie et l'économie. La nécessité de résoudre ces équations, souvent sans solutions analytiques, a donné naissance à des approches numériques sophistiquées. Face à la complexité croissante des systèmes à modéliser, il est essentiel d'acquérir des compétences techniques tout en comprenant l'impact des méthodes numériques sur la résolution de problèmes réels. Ce cours vise à établir un lien entre la théorie mathématique et son application pratique, préparant ainsi à des défis futurs dans la carrière.

Le projet vise à modéliser et résoudre numériquement des problèmes d'advection-diffusion-réaction pour un scalaire passif en mettant en œuvre des techniques avancées de discrétisation et d'adaptation de maillage. Les objectifs incluent :

### Modélisation : 
Développement d’un modèle d’advection-diffusion-réaction pour des solutions stationnaires et instationnaires.
### Estimation d’Erreur :
Utilisation de méthodes a posteriori pour estimer l'erreur de discrétisation.
### Adaptation de Maillage : 
Application de contrôles de métrique pour adapter dynamiquement le maillage en fonction de l'erreur locale.
### Optimisation de Contrôle : 
Optimisation de la distribution de sources pour un scalaire passif en contrôlant les sources pour minimiser l'erreur par rapport à une solution cible.

---

## Séance 1 : Introduction à l'Environnement de Développement

### Objectif
Familiariser avec l'environnement de développement Python, préparant ainsi le terrain pour la résolution des EDP.

### Activités Détaillées
- **Installation et Configuration** : 
  - **Tâche** : Installer Spyder et Jupyter Notebook sur chaque machine.
  - **Développement** : Suivre les instructions pas à pas pour garantir le bon fonctionnement de tous les logiciels, incluant le téléchargement des bibliothèques nécessaires, telles que NumPy et Matplotlib.

- **Exploration de l'Interface** :
  - **Tâche** : Se familiariser avec l'interface de Spyder et Jupyter.
  - **Développement** : Explorer les différentes fonctionnalités, apprendre à exécuter des scripts, à utiliser le terminal et à visualiser les résultats, ce qui permet une navigation efficace durant les travaux pratiques.

- **Introduction à l'Équation à Résoudre** : 
  - **Tâche** : Présenter l'équation générale.
  - **Développement** : Décortiquer l'équation \( u_t + \langle V, \nabla u \rangle - \nu \Delta u = -\lambda u + f \). Apprendre les termes impliqués, tels que les conditions aux limites et les implications physiques de chaque composante.

- **Tâches à Réaliser** :
  - **Tâche** : Utiliser l'IA générative pour créer un code Python pour résoudre l'équation.
  - **Développement** : Développer un code de simulation en définissant des maillages et les intégrant dans les algorithmes. En résolvant une version simplifiée \( u_t = -\lambda u \), comparer les résultats numériques avec la solution analytique.

---

## Séance 2 : Analyse du Code `adrs.py`

### Contexte
Cette séance se concentre sur l'analyse du code `adrs.py`, un outil essentiel pour résoudre des EDP à une dimension.

### Objectifs
1. Développer une compréhension approfondie des composants du code.
2. Identifier et déterminer la fonction \( f(s) \) pour les simulations.

### Activités Détaillées
- **Lecture du Code `adrs.py`** :
  - **Tâche** : Analyser les différents blocs du code.
  - **Développement** : Identifier les algorithmes utilisés pour résoudre l'équation, en portant une attention particulière à la structure du code, aux boucles d'itération et à la gestion des conditions aux limites.

- **Équation à Résoudre** : 
  - **Tâche** : Formuler l'équation avec les conditions appropriées.
  - **Développement** : Formuler l'équation avec des conditions initiales, modélisant des situations réelles par le biais de simulations numériques.

- **Implémentation de \( f(s) \)** :
  - **Tâche** : Déterminer la fonction \( f(s) \) pour la simulation.
  - **Développement** : La définition précise de \( f(s) \) est cruciale pour l'exactitude des résultats, en discutant des implications de différentes fonctions.

- **Vérification de la Convergence** :
  - **Tâche** : Tester le code avec un maillage de 100 points.
  - **Développement** : Analyser les résultats obtenus et discuter de la convergence de la solution, en apprenant à ajuster les paramètres pour améliorer les résultats.

- **Visualisation** :
  - **Tâche** : Tracer la solution numérique.
  - **Développement** : Comparer visuellement la solution numérique avec la solution exacte, en calculant l'erreur L2 pour évaluer la performance des algorithmes.

---

## Séances 3-4 : Adaptation A Posteriori

### Thématique
Exploration des méthodes d'adaptation a posteriori pour optimiser la précision des solutions numériques.

### Objectifs
1. Appliquer l'adaptation a posteriori dans le cadre de l'intégration.
2. Comparer l'efficacité de l'intégration de Riemann et de Lebesgue.

### Activités Détaillées
- **Analyse de la Fonction \( f(x) \)** :
  - **Tâche** : Présenter et analyser la fonction à intégrer.
  - **Développement** : Examiner la fonction \( f(x) = \exp(-a(x-L/3)^2) - 2\exp(-a(x-2L/3)^2) \). Discuter de son comportement et de son importance dans les méthodes d'intégration.

- **Intégration de Riemann** :
  - **Tâche** : Implémenter une méthode d'intégration.
  - **Développement** : Utiliser des pas non uniformes et surmonter les défis liés à la gestion des discontinuités, ouvrant la discussion sur la robustesse des méthodes d'intégration.

- **Adaptation par Lebesgue** :
  - **Tâche** : Appliquer l'intégration de Lebesgue.
  - **Développement** : Mettre en œuvre cette méthode pour un pas uniforme et évaluer l'erreur d'intégration, apprenant à comparer les performances des deux méthodes.

- **Comparaison des Méthodes** :
  - **Tâche** : Analyser les résultats des deux méthodes.
  - **Développement** : Une discussion enrichissante sur les avantages et inconvénients de chaque méthode en termes de précision et de complexité.

---

## Séance 5 : Résolution Instationnaire

### Contexte
Utilisation du code `adrs_insta.py` pour traiter des équations instationnaires, avec une attention particulière à la visualisation des résultats.

### Objectifs
1. Analyser les erreurs à différents instants pour des maillages uniformes.
2. Établir des critères d'arrêt pour les itérations d'adaptation.

### Activités Détaillées
- **Analyse du Code `adrs_insta.py`** :
  - **Tâche** : Étudier la structure du code pour les EDP instationnaires.
  - **Développement** : Comprendre les itérations temporelles, les conditions aux limites et la mise en œuvre de l'adaptation a posteriori dans le code.

- **Simulation de l'Équation** :
  - **Tâche** : Résoudre l'équation \( u_t + \Delta u = 0 \) à l'aide de `adrs_insta.py`.
  - **Développement** : Effectuer des simulations pour différentes valeurs de paramètres, en analysant l'impact sur la stabilité et la précision des solutions.

- **Visualisation des Résultats** :
  - **Tâche** : Tracer les solutions à différents instants.
  - **Développement** : Créer des graphiques illustrant l'évolution de la solution, en mettant en avant les différences entre les solutions numériques et analytiques.

- **Analyse des Erreurs** :
  - **Tâche** : Calculer les erreurs L2 pour différents maillages.
  - **Développement** : Discuter des résultats et des implications de la résolution des erreurs pour les méthodes d'adaptation.

---

## Séance 6 : Optimisation et Visualisation des Résultats

### Contexte
Cette séance se concentre sur l'optimisation des paramètres d'algorithme et sur la visualisation avancée des résultats numériques.

### Objectifs
1. Améliorer la précision des simulations par l'optimisation des paramètres.
2. Développer des techniques de visualisation pour mieux interpréter les résultats.

### Activités Détaillées
- **Optimisation des Paramètres** :
  - **Tâche** : Identifier les paramètres clés à optimiser dans les simulations.
  - **Développement** : Expérimenter avec différentes configurations de maillage et de pas temporels, en observant l'impact sur la stabilité et la convergence des résultats.

- **Techniques de Visualisation** :
  - **Tâche** : Appliquer des méthodes de visualisation avancées, telles que les heatmaps et les surfaces 3D.
  - **Développement** : Utiliser des bibliothèques comme Matplotlib et Seaborn pour créer des visualisations interactives, permettant une meilleure interprétation des données.

- **Discussion sur les Résultats** :
  - **Tâche** : Analyser et discuter des résultats obtenus.
  - **Développement** : Évaluer l'impact des choix d'optimisation et des techniques de visualisation sur la compréhension des phénomènes modélisés.

---
## Résultats
Les scripts de ce dépôt génèrent divers graphiques et visualisations, notamment :

### Convergence :
Analyse de la convergence des solutions numériques par rapport aux solutions exactes.
### Évolution de l'erreur : 
Visualisation de l'évolution de l'erreur L2 et 𝐻1 en fonction du pas de maillage.
### Adaptation de Maillage : 
Comparaison de la distribution du maillage pour les solutions stationnaires et instationnaires à différents instants.
### Optimisation de Contrôle :
Visualisation des résultats de l’optimisation de la distribution des sources pour atteindre une solution cible.


---
## Conclusion
Au terme de ces séances, une compréhension approfondie des méthodes d'analyse numérique pour les EDP a été développée. Les compétences acquises dans l'utilisation des outils de programmation, d'adaptation de maillage et d'optimisation des paramètres constituent un atout essentiel pour résoudre des problèmes complexes. La capacité à allier théorie et pratique renforce la préparation à des défis futurs dans le domaine des mathématiques appliquées.

Ces séances ont non seulement permis d'acquérir des compétences techniques, mais ont également favorisé un esprit critique face aux méthodes utilisées. La discussion autour des choix méthodologiques, des erreurs observées et des résultats numériques a enrichi l'expérience d'apprentissage. En outre, l'accent mis sur la visualisation des résultats a renforcé l'importance de la communication des résultats scientifiques, une compétence cruciale dans le domaine de la recherche.

L'interaction entre théorie et pratique, ainsi que l'application des concepts mathématiques à des problèmes concrets, a préparé les participants à aborder des travaux futurs avec une confiance accrue. Ces compétences seront sans aucun doute précieuses dans la poursuite d'une carrière dans les mathématiques appliquées ou dans des disciplines connexes. 

