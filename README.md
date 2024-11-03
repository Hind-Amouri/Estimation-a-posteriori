# Estimation-a-posteriori

# Motivations
Ce projet est basé sur un cours qui se concentre sur la modélisation numérique à l'aide des équations aux dérivées partielles (EDP).Les principaux objectifs sont les suivants :

Modélisation d’un problème d’advection-diffusion-réaction-production pour un scalaire passif.
Discrétisation du problème sur un maillage à pas variable pour des solutions stationnaires.
Estimation de l’erreur de discrétisation.
Adaptation du maillage en utilisant un contrôle de métrique qui ajuste la distorsion locale de l’espace.
Extension de ces concepts aux problèmes instationnaires.
Optimisation de la distribution du scalaire passif en contrôlant ses sources.

L'objectif final est de mettre en œuvre ces concepts dans des programmes Python et d'explorer les méthodes d'adaptation de maillage pour améliorer la précision des solutions numériques.

# Séances
# Séance 1 : Prise en Main de Spyder/Jupyter et Résolution Numérique
Cette première séance est dédiée à la prise en main de Spyder ou Jupyter Notebook pour les calculs scientifiques. L’objectif principal est la résolution d’une équation d’advection-diffusion dans un domaine rectangulaire, avec un schéma numérique simple.

Équation modélisée :
𝑢,𝑡+𝑣1𝑢,𝑥−𝑛𝑢⋅𝑢,𝑥𝑥=−𝜆𝑢+𝑓u,t+v1u,x−nu⋅u,xx=−λu+f

But : Implémenter un schéma d’Euler explicite pour une équation différentielle simple :
𝑢,𝑡=−𝜆𝑢,𝑢(0)=1,𝜆=1u,t=−λu,u(0)=1,λ=1

Résultat attendu : Comparaison des solutions numériques et exactes pour différentes valeurs du pas de temps.

# Séance 2 : Analyse de la Convergence et Adaptation de Maillage
Durant cette séance, nous analysons la convergence de la solution numérique en fonction du maillage et mettons en place un critère d’adaptation de maillage pour améliorer la précision.

Problème traité : Résolution d’un problème d’advection-diffusion avec 
𝑢(𝑠)=exp⁡(−10(𝑠−𝐿/2)2)u(s)=exp(−10(s−L/2) 2 ), et calcul des erreurs 
𝐿2L2 et 𝐻1H1.
Objectif : Adapter le maillage pour obtenir la meilleure solution possible en fonction de l'erreur de discrétisation.

# Séance 3-4 : Adaptation de Maillage à Posteriori
Nous comparons deux méthodes d'intégration (Riemann et Lebesgue) et introduisons un maillage adaptatif en fonction de l'erreur.

Fonction traitée :
𝑓(𝑥)=exp⁡(−𝑎(𝑥−𝐿/3)2)−2exp(−𝑎(𝑥−2𝐿/3)2),𝑎=100,𝐿=1f(x)=exp(−a(x−L/3) 2 )−2exp(−a(x−2L/3) 2 ),a=100,L=1

Objectif : Implémenter un maillage adaptatif basé sur la fonction 𝑓(𝑥)
f(x), et tracer l’évolution du nombre de points de maillage en fonction de l’erreur de discrétisation.
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

# **Séance 6 : Optimisation et Contrôle Numérique**

Dans cette séance, nous explorons l'optimisation de la distribution des sources dans un problème d'advection-diffusion-réaction en utilisant un modèle linéaire pour simplifier les calculs.

## **Code `optim_adrs.py` : Optimisation de la Distribution des Sources**
- **Objectif** : Déterminer une distribution optimale de sources pour obtenir une solution ciblée, \( u_{\text{des}} \), en minimisant l'erreur.
- **Description du problème** :
  - Formulation du problème inverse en définissant un objectif \( J(x) = \frac{1}{2} ||u(x) - u_{\text{des}}||^2 \).
  - Utilisation de la linéarité de l'équation ADRS pour simplifier les calculs des gradients et de la matrice Hessienne.

### **Algorithme**
1. Calculer une solution cible \( u_{\text{des}} \) pour un vecteur de contrôle donné \( x_{\text{cible}} \).
2. Minimiser \( J(x) \) en résolvant \( A x = b \) avec \( A_{ij} = \langle u_i, u_j \rangle_{L2} \) et \( b_i = \langle u_i, u_{\text{des}} - u_0 \rangle_{L2} \).
3. Utiliser un algorithme de minimisation pour obtenir le vecteur de contrôle optimal \( x_{\text{opt}} \).
4. Introduire une boucle de raffinement du maillage pour améliorer la convergence de la solution vers \( x_{\text{opt}} \).

### **Visualisations et Résultats**
- **Comparaisons de Contrôle** : Visualisation de la surface \( J(x_1, x_2) \) en échantillonnant les deux premiers contrôles, les autres étant fixes.
- **Impact de l'Adaptation de Maillage** : Comparaison de l'erreur obtenue avec et sans adaptation de maillage.
- **Convergence du Contrôle** : Analyse de la convergence de \( x_{\text{opt}} \) avec le raffinement du maillage.

## Prérequis
- Python 3.x
- `numpy`, `matplotlib`, `scipy`

## Exécution
Pour exécuter le script :
```bash
python optim_adrs.py

## **Structure du Projet**
Le projet est organisé autour de plusieurs fichiers Python qui implémentent les concepts abordés dans les séances :

- `adrs.py` : Code de base pour la résolution de l'équation d’advection-diffusion.
- `adrs_multiple_mesh_adap.py` : Adaptation de maillage pour des solutions stationnaires.
- `adrs_insta_multiple_mesh_adap.py` : Adaptation instationnaire avec un terme source dépendant du temps.

## **Installation**
Pour exécuter les codes fournis, vous devez installer les bibliothèques suivantes :

- Python 3.x
- `numpy` pour les calculs numériques
- `matplotlib` pour la visualisation
- Jupyter Notebook ou Spyder pour l’exécution interactive
