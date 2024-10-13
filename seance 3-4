## Séances 3 et 4 : Adaptation à Posteriori et Intégration

### Objectifs :
- Comprendre l'adaptation à posteriori en comparant les méthodes d'intégration de Riemann et Lebesgue.
- Implémenter l'intégration de la fonction \( f(x) = \exp(-a(x - L/3)^2) - 2\exp(-a(x - 2L/3)^2) \) avec \( a = 100 \) et \( L = 1 \).
- Mettre en place un maillage adaptatif et un critère d'arrêt mixte (nombre de points et erreur \( L2 \)) dans le code `adrs_multiple_mesh_adap.py`.

### Méthodologie :
1. **Intégration de Riemann et Lebesgue :**
   - Nous avons comparé deux méthodes d'intégration pour calculer l'intégrale de la fonction donnée.
   - La méthode de Riemann utilise un pas uniforme en \( x \), tandis que la méthode de Lebesgue utilise un pas uniforme en \( y = f(x) \).

2. **Adaptation du maillage :**
   - Le code `adrs_multiple_mesh_adap.py` a été restructuré pour introduire trois fonctions principales : `adrs()`, `metric()`, et `mesh()`.
   - Nous avons mis en place un critère d'arrêt mixte basé sur le nombre de points de maillage et l'erreur \( L2 \).

3. **Résultats :**
   - Les courbes de \( NX(\text{err}) \) pour différentes valeurs d'erreur montrent comment le nombre de points de maillage évolue en fonction de l'erreur tolérée.
