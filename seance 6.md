# Séance 6 : Optimisation et Contrôle Numérique

Dans cette séance, nous explorons l'optimisation de la distribution des sources dans un problème d'advection-diffusion-réaction en utilisant un modèle linéaire pour simplifier les calculs.

## Contenu

### **Code `optim_adrs.py` : Optimisation de la Distribution des Sources**
- **Objectif** : Déterminer une distribution optimale de sources pour obtenir une solution ciblée, \( u_{\text{des}} \), en minimisant l'erreur.
- **Description du problème** :
  - Formulation du problème inverse en définissant un objectif \( J(x) = \frac{1}{2} ||u(x) - u_{\text{des}}||^2 \).
  - Utilisation de la linéarité de l'équation ADRS pour simplifier les calculs des gradients et de la matrice Hessienne.
  
- **Algorithme** :
  1. Calculer une solution cible \( u_{\text{des}} \) pour un vecteur de contrôle donné \( x_{\text{cible}} \).
  2. Minimiser \( J(x) \) en résolvant \( A x = b \) avec \( A_{ij} = \langle u_i, u_j \rangle_{L2} \) et \( b_i = \langle u_i, u_{\text{des}} - u_0 \rangle_{L2} \).
  3. Utiliser un algorithme de minimisation pour obtenir le vecteur de contrôle optimal \( x_{\text{opt}} \) en fonction de la solution cible.
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
