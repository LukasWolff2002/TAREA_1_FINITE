
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def plot_original_structure_all_forces(nodes, elements):
    fig, ax = plt.subplots(figsize=(8, 6))

    # Obtener todas las áreas para normalización de colores
    areas = np.array([element.A for element in elements])
    norm = mcolors.Normalize(vmin=np.min(areas), vmax=np.max(areas))
    cmap = cm.viridis

    # Graficar elementos originales con color basado en el área
    for element in elements:
        x1, y1 = element.n1.coord
        x2, y2 = element.n2.coord
        color = cmap(norm(element.A))
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=2, label="Elemento (color por A)" if element == elements[0] else "")

    # Graficar nodos
    for node in nodes:
        x, y = node.coord
        ax.scatter(x, y, color='g', s=50, zorder=3, label="Nodo" if node == nodes[0] else "")

    # Graficar todas las fuerzas aplicadas en cada nodo (Fuerzas en X, Y)
    for node in nodes:
        x, y = node.coord
        fx, fy, m = node.force_vector  # Considerar también el momento

        if fx != 0:
            ax.quiver(x, y, 1, 0, angles='xy', scale_units='xy', scale=1, color='b',
                      label="Fuerza en X" if node == nodes[0] else "")
            ax.text(x + 0.3, y, f'Fx={fx:.2f}', fontsize=9, color='b')

        if fy != 0:
            ax.quiver(x, y, 0, -1, angles='xy', scale_units='xy', scale=1, color='r',
                      label="Fuerza en Y" if node == nodes[0] else "")
            ax.text(x, y - 0.3, f'Fy={fy:.2f}', fontsize=9, color='r')

        if m != 0:
            arc = plt.Circle((x, y), 0.2, color='purple', fill=False, linestyle='dashed')
            ax.add_patch(arc)
            ax.text(x + 0.25, y + 0.25, f'M={m:.2f}', fontsize=10, color='purple')

    # Barra de color para indicar valores de área
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Configuración del gráfico
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title("Estructura Original con Colores según Área de la Sección")
    ax.axis("equal")
    plt.grid(True)
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def plot_deformed_structure(nodes, elements, scale=1000):
    """
    Grafica la estructura deformada usando funciones de forma de viga Euler-Bernoulli 
    y transformando correctamente los desplazamientos locales a coordenadas globales,
    según el ángulo del elemento.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Graficar estructura original
    for element in elements:
        x1, y1 = element.n1.coord
        x2, y2 = element.n2.coord
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1, label="Original" if element == elements[0] else "")

    # Graficar estructura deformada con funciones de forma
    for element in elements:
        xi = element.n1.coord
        xf = element.n2.coord
        L = element.L
        angle = element.angle
        c = np.cos(angle)
        s = np.sin(angle)

        # Desplazamientos nodales en local: [u1, v1, θ1], [u2, v2, θ2]
        u1, v1, th1 = element.n1.def_vector
        u2, v2, th2 = element.n2.def_vector

        # Número de puntos a interpolar
        n_points = 40
        x_local = np.linspace(0, L, n_points)

        # Funciones de forma para flexión (Euler-Bernoulli)
        def shape_functions(x, L):
            ξ = x / L
            N1 = 1 - 3*ξ**2 + 2*ξ**3
            N2 = x * (1 - 2*ξ + ξ**2)
            N3 = 3*ξ**2 - 2*ξ**3
            N4 = x * (ξ - 1)**2
            return N1, N2, N3, N4

        # Matriz de rotación local → global
        R = np.array([[c, -s],
                      [s,  c]])

        # Lista de puntos transformados
        xg_list, yg_list = [], []

        for x in x_local:
            N1, N2, N3, N4 = shape_functions(x, L)
            v = N1 * v1 + N2 * th1 + N3 * v2 + N4 * th2  # desplazamiento vertical local (v(x))
            local_point = np.array([x, v * scale])      # deformación amplificada
            global_point = xi + R @ local_point         # transformación completa
            xg_list.append(global_point[0])
            yg_list.append(global_point[1])

        ax.plot(xg_list, yg_list, 'r--', linewidth=1.5, label="Deformada" if element == elements[0] else "")

    # Graficar nodos originales y desplazados
    for node in nodes:
        x, y = node.coord
        dx, dy = node.def_vector[:2] * scale
        ax.scatter(x, y, color='g', s=30, zorder=3, label="Nodo Original" if node == nodes[0] else "")
        ax.scatter(x + dx, y + dy, color='b', s=30, zorder=3, label="Nodo Deformado" if node == nodes[0] else "")

    ax.set_title("Estructura Original y Deformada (con ángulo y funciones de forma)")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.axis('equal')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()
