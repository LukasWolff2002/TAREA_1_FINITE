
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

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

def plot_structure_with_local_displacements(nodes, elements, scale=1000):
    """
    Graficar la estructura usando los desplazamientos locales de cada elemento en el sistema local.
    Los desplazamientos locales de cada nodo se transforman a coordenadas globales y se grafican.
    
    Parámetros:
    - nodes: Lista de nodos de la estructura.
    - elements: Lista de elementos (con .angle, .L, .n1, .n2 y .ug).
    - scale: Factor para amplificar la deformación visualmente.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Graficar estructura original
    for element in elements:
        x1, y1 = element.n1.coord
        x2, y2 = element.n2.coord
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1, label="Original" if element == elements[0] else "")

    # Graficar deformada usando desplazamientos locales
    for element in elements:
        xi = element.n1.coord
        xf = element.n2.coord
        L = element.L
        angle = element.angle
        c = np.cos(angle)
        s = np.sin(angle)

        # Desplazamientos locales de los nodos [uL, vL, θL] de cada nodo
        u1_L, v1_L, th1_L = element.ug[0], element.ug[1], element.ug[2]  # Nodo 1
        u2_L, v2_L, th2_L = element.ug[3], element.ug[4], element.ug[5]  # Nodo 2

        # Función para transformar desplazamientos locales a globales
        def transform_to_global(u_L, v_L, th_L, angle):
            R = np.array([[c, s, 0], 
                          [-s, c, 0], 
                          [0, 0, 1]])  # Matriz de rotación 2D

            local_displacement = np.array([u_L, v_L, th_L])
            global_displacement = R @ local_displacement
            return global_displacement

        # Transformación de desplazamientos locales a globales
        u1_G, v1_G, th1_G = transform_to_global(u1_L, v1_L, th1_L, angle)
        u2_G, v2_G, th2_G = transform_to_global(u2_L, v2_L, th2_L, angle)

        # Número de puntos a interpolar
        n_points = 40
        x_local = np.linspace(0, L, n_points)

        # Elemento horizontal (viga): usamos las funciones de forma de viga de Euler-Bernoulli
        if angle != np.pi / 2 and angle != 0:  # No es un elemento vertical
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
                v = N1 * v1_G + N2 * th1_G + N3 * v2_G + N4 * th2_G  # desplazamiento vertical local (v(x))
                v_scaled = v * scale

                local_point = np.array([x, v_scaled])      # deformación amplificada
                global_point = xi + R @ local_point         # transformación a sistema global
                xg_list.append(global_point[0])
                yg_list.append(global_point[1])

            ax.plot(xg_list, yg_list, 'r--', linewidth=1.5, label="Deformada" if element == elements[0] else "")

        else:
            # Elemento vertical: rotar 90 grados y considerar deformación solo en el eje vertical
            dx = xi[0]  # X constante para elementos verticales
            dy = np.linspace(xi[1], xf[1], n_points)  # Interpolamos el desplazamiento vertical

            # Rotación de 90 grados: intercambiar las coordenadas X y Y
            x_rot = xi[1]
            ax.plot(np.full_like(dy, x_rot), dy, 'r--', linewidth=1.5, label="Deformada Vertical" if element == elements[0] else "")

    # Graficar nodos originales y desplazados
    for node in nodes:
        x, y = node.coord
        dx, dy = node.def_vector[:2] * scale
        ax.scatter(x, y, color='g', s=40, zorder=3, label="Nodo Original" if node == nodes[0] else "")
        ax.scatter(x + dx, y + dy, color='b', s=40, zorder=3, label="Nodo Deformado" if node == nodes[0] else "")

    ax.set_title("Estructura con Desplazamientos Locales (Deformada vs Original)")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.axis('equal')
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    plt.show()


def plot_deformed_structure(nodes, elements, scale=1000):
    """
    Graficar la estructura deformada usando las funciones de forma para cada elemento,
    combinando los desplazamientos axiales y de flexión. La deformación se calcula
    en función de las posiciones de los nodos de cada elemento.
    
    Parámetros:
    - nodes: Lista de nodos de la estructura.
    - elements: Lista de elementos (con .angle, .L, .n1, .n2 y .ug).
    - scale: Factor para amplificar la deformación visualmente.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    # Graficar la estructura original
    for element in elements:
        x1, y1 = element.n1.coord
        x2, y2 = element.n2.coord
        #ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1)  # Original structure

    # Graficar la estructura deformada
    for element in elements:
        ug = element.ug  # Desplazamientos locales

        # Condiciones de borde para comportamiento axial
        u0 = 0  # Desplazamiento axial inicial en el nodo 1
        u1 = ug[3]  # Desplazamiento axial en el nodo 2

        L = element.L  # Longitud del elemento
        x_local = np.linspace(0, L, 1000)  # Puntos a lo largo de la longitud del elemento

        # Funciones de forma para el comportamiento axial
        n0 = 1 - (x_local / L)
        n1 = x_local / L

        # Desplazamientos axiales a lo largo del elemento
        u_axial = n0 * u0 + n1 * u1  # Desplazamiento axial en cada punto

        # Funciones de forma para la flexión (viga de Euler-Bernoulli)
        u2 = ug[1]  # Desplazamiento vertical del nodo 1
        u3 = ug[2]  # Rotación del nodo 1
        u4 = ug[4]  # Desplazamiento vertical del nodo 2
        u5 = ug[5]  # Rotación del nodo 2

        # Funciones de forma para flexión
        def shape_functions_flexion(x, L):
            ξ = (2 * x / L) - 1  # Función de forma
            N1 = (1 / 4) * ((1 - ξ) ** 2) * (2 + ξ) * u2
            N2 = (L / 8) * ((1 - ξ) ** 2) * (1 + ξ) * u3
            N3 = (1 / 4) * ((1 + ξ) ** 2) * (2 - ξ) * u4
            N4 = -(L / 8) * ((1 + ξ) ** 2) * (1 - ξ) * u5
            return N1, N2, N3, N4

        # Desplazamientos a lo largo del elemento (de flexión)
        v_flexion = np.zeros_like(x_local)

        for i, x in enumerate(x_local):
            N1, N2, N3, N4 = shape_functions_flexion(x, L)
            v_flexion[i] = N1 + N2 + N3 + N4  # Desplazamiento vertical total debido a flexión

        # Desplazamiento total (combinando axial y flexión)
        u_total = u_axial + v_flexion * scale  # Deformación total combinada (axial + flexión)

        # Calcular las posiciones finales de los puntos de la estructura deformada
        xi, yi = element.n1.coord
        xf, yf = element.n2.coord

        # Graficar la deformación combinada (axial + flexión) para este elemento
        # Asegurándonos de que los desplazamientos estén relativos a los nodos
        if xi == xf:  # Elemento vertical
            print(yi)
            # Si el elemento es vertical, tratamos la deformación como si fuera horizontal
            # Primero, calculamos la deformación como si fuera un elemento horizontal
            x_deformed = xi + np.zeros_like(x_local)  # La coordenada X permanece constante
            y_deformed = yi + u_total  # Deformación en el eje Y debido a la flexión y el desplazamiento axial
            ax.plot(x_deformed, y_deformed, 'r--', linewidth=1.5)
        else:
            # Elemento horizontal (conforme a la deformación en X e Y)
            x_deformed = xi + (xf - xi) * x_local / L  # Ajuste de la posición en X (relativo a nodos)
            y_deformed = yi + u_total  # Ajuste de la posición en Y (deformación relativa)
            ax.plot(x_deformed, y_deformed, 'r--', linewidth=1.5)

    ax.set_title("Estructura Deformada (Axial + Flexión)")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.axis('equal')
    ax.grid(True)
    plt.tight_layout()
    plt.show()