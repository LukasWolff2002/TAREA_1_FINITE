import numpy as np
import matplotlib.pyplot as plt

# Parámetros de la viga
L = 10  # Longitud de la viga (m)

# Definir el cambio de variable chi
chi_vals = np.linspace(-1, 1, 100)

# Funciones de forma
f0 = (1/4) * (1 - chi_vals)**2 * (2 + chi_vals)
f1 = (L / 8) * (1 - chi_vals)**2 * (1 + chi_vals)
f2 = (1/4) * (1 + chi_vals)**2 * (2 - chi_vals)
f3 = -(L / 8) * (1 + chi_vals)**2 * (1 - chi_vals)

# Crear una figura y una cuadrícula de subgráficos
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Graficar la primera función de forma (f0)
axs[0, 0].plot(chi_vals, f0, color='b')
axs[0, 0].set_title(r'$f_0(\chi) = \frac{1}{4}(1-\chi)^2(2+\chi)$')
axs[0, 0].set_xlabel(r'$\chi$')
axs[0, 0].set_ylabel('Función de forma')
axs[0, 0].grid(True)

# Graficar la segunda función de forma (f1)
axs[0, 1].plot(chi_vals, f1, color='g')
axs[0, 1].set_title(r'$f_1(\chi) = \frac{L}{8}(1-\chi)^2(1+\chi)$')
axs[0, 1].set_xlabel(r'$\chi$')
axs[0, 1].set_ylabel('Función de forma')
axs[0, 1].grid(True)

# Graficar la tercera función de forma (f2)
axs[1, 0].plot(chi_vals, f2, color='r')
axs[1, 0].set_title(r'$f_2(\chi) = \frac{1}{4}(1+\chi)^2(2-\chi)$')
axs[1, 0].set_xlabel(r'$\chi$')
axs[1, 0].set_ylabel('Función de forma')
axs[1, 0].grid(True)

# Graficar la cuarta función de forma (f3)
axs[1, 1].plot(chi_vals, f3, color='purple')
axs[1, 1].set_title(r'$f_3(\chi) = -\frac{L}{8}(1+\chi)^2(1-\chi)$')
axs[1, 1].set_xlabel(r'$\chi$')
axs[1, 1].set_ylabel('Función de forma')
axs[1, 1].grid(True)

# Ajustar el espacio entre los subgráficos
plt.tight_layout()

# Mostrar los gráficos
plt.show()
