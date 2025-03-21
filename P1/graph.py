
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


def plot_deformed_structure(nodes, elements, scale=1000):
    """
    Grafica la estructura desplazada usando los valores de desplazamiento almacenados en cada nodo.
    
    Parámetros:
    - nodes: Lista de nodos de la estructura.
    - elements: Lista de elementos que conectan los nodos.
    - scale: Factor de escala para amplificar las deformaciones y hacerlas visibles.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Diccionarios para almacenar coordenadas originales y desplazadas
    coords_original = {node.n: node.coord for node in nodes}
    coords_deformed = {}

    # Calcular coordenadas desplazadas
    for node in nodes:
        dx, dy = node.def_vector[:2] * scale  # Aplicar factor de escala
        coords_deformed[node.n] = node.coord + np.array([dx, dy])

    # Graficar estructura original
    for element in elements:
        x1, y1 = coords_original[element.n1.n]
        x2, y2 = coords_original[element.n2.n]
        ax.plot([x1, x2], [y1, y2], 'k-', label="Original" if element == elements[0] else "")

    # Graficar estructura desplazada
    for element in elements:
        xd1, yd1 = coords_deformed[element.n1.n]
        xd2, yd2 = coords_deformed[element.n2.n]
        ax.plot([xd1, xd2], [yd1, yd2], 'r--', label="Deformada" if element == elements[0] else "")

    # Graficar nodos originales y desplazados
    for node in nodes:
        x, y = coords_original[node.n]
        ax.scatter(x, y, color='g', s=50, zorder=3, label="Nodo Original" if node == nodes[0] else "")

        xd, yd = coords_deformed[node.n]
        ax.scatter(xd, yd, color='b', s=50, zorder=3, label="Nodo Desplazado" if node == nodes[0] else "")

    # Configuración del gráfico
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title("Estructura Original y Deformada")
    ax.legend()
    ax.axis("equal")
    plt.grid(True)
    plt.show()

