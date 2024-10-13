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

# Séance 4-5 : Adaptation Instationnaire avec Terme Source Temporel
Nous introduisons un terme source dépendant du temps dans le problème d'advection-diffusion, ce qui conduit à une solution instationnaire.

Équation modélisée :
𝑢(𝑡)𝑣(𝑠)=sin⁡(4𝜋𝑡)𝑣(𝑠)u(t)v(s)=sin(4πt)v(s)
Objectif : Adapter le maillage en temps et en espace pour tenir compte des variations temporelles de la solution, tout en contrôlant l’erreur et en optimisant le nombre de points de maillage.

# Structure du Projet
Le projet est organisé autour de plusieurs fichiers Python qui implémentent les concepts abordés dans les séances :

adrs.py : Code de base pour la résolution de l'équation d’advection-diffusion.

adrs_multiple_mesh_adap.py : Adaptation de maillage pour des solutions stationnaires.

adrs_insta_multiple_mesh_adap.py : Adaptation instationnaire avec un terme source dépendant du temps.

Installation
Pour exécuter les codes fournis, vous devez installer les bibliothèques suivantes :

Python 3.x
Numpy pour les calculs numériques
Matplotlib pour la visualisation
Jupyter Notebook ou Spyder pour l’exécution interactive
Installez les dépendances avec la commande suivante :


Résultats
Le projet génère des graphiques et des visualisations montrant :

La convergence des solutions numériques vers les solutions exactes.
L'évolution de l'erreur 
L2 et H1 en fonction du pas de maillage.
Des visualisations des solutions instationnaires à différents instants de temps.
