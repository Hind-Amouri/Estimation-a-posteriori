import numpy as np
import matplotlib.pyplot as plt

# Paramètres
lambda_ = 1
u0 = 1
Dt = 1  # pas de temps
T = 60  # durée en secondes
N = int(T/Dt)  # nombre de points de temps

# Initialisation
t = np.linspace(0, T, N+1)
u_exact = u0 * np.exp(-lambda_ * t)
u_numeric = np.zeros(N+1)
u_numeric[0] = u0

# Schéma d'Euler explicite
for n in range(N):
    u_numeric[n+1] = u_numeric[n] - Dt * lambda_ * u_numeric[n]

# Calcul de l'erreur L2
error = np.sqrt(np.sum((u_exact - u_numeric)**2) * Dt)

# Tracer les solutions
plt.figure(figsize=(10, 6))
plt.plot(t, u_exact, label="Solution exacte", linestyle="--", color="blue")
plt.plot(t, u_numeric, label="Solution numérique", color="red")
plt.xlabel('Temps (s)')
plt.ylabel('u(t)')
plt.legend()
plt.title('Comparaison des solutions exactes et numériques')
plt.grid(True)
plt.show()

# Tracer l'erreur
plt.figure(figsize=(10, 6))
plt.plot(t, u_exact - u_numeric, label="Erreur", color="green")
plt.xlabel('Temps (s)')
plt.ylabel('Erreur')
plt.title('Erreur en fonction du temps')
plt.grid(True)
plt.show()

print(f"L'erreur L2 de la solution numérique est : {error}")
