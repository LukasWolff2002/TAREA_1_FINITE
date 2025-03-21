import numpy as np
from nodes import Node
from elements import Elements
from desplazamientos import Desplazamientos
from assembly import Assembly
import matplotlib.pyplot as plt

nodes = []
elements = []

espaciado_h = 7 #m
Espaciado_v = [3.66, 
               5.49, 
               (3.96)/2, 
               (3.96)/2, 
               3.96, 
               (3.96)/2, 
               (3.96/2), 
               3.96, 
               (3.96)/2, 
               (3.96/2), 
               3.96, 
               (3.96)/2, 
               (3.96/2),
               3.96]

#Defino los pisos en que hay cambios de seccion
losas = [1, #Piso 1
         1, #Piso 2
         0, #Piso 3
         1, #Piso 4
         1, #Piso 5
         0, #Piso 6
         1, #Piso 7
         1, #Piso 8
         0, #Piso 9
         1, #Piso 10
         1, #Piso 11
         0, #Piso 12
         1, #Piso 13
        1] #Piso 14
        
#Ahora defino las secciones por piso

        #Seccion   Izquierda   Centrales    Derecha
secciones_piso = [[[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 1
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 2
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 3
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 4
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 5
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 6
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 7
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 8
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 9
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 10
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 11
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 12
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]], #Piso 13
                  [[0.3, 0.3], [1.1, 1.1], [0.1, 0.1]]  #Piso 14
                  ]

nodos_ancho = 6
pisos = 14

base_v = 0
#Defino los nodos base
for i in range(nodos_ancho):
    nodes.append(Node(i, np.array([i*espaciado_h, 0]), np.array([1, 1, 1]), np.array([0.0, 0.0, 0.0])))


#Defino los nodos de los pisos
for j in range(1, pisos+1):
    espaciado_v = Espaciado_v[j-1]
    vertical = base_v + espaciado_v

   
    for i in range(nodos_ancho):
        nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

    #Ahora conecto verticalmente los pisos
    for i in range(nodos_ancho):
        
        if i == 0:
            A = secciones_piso[j-1][0]
        elif i == nodos_ancho-1:
            A = secciones_piso[j-1][2]
        else:
            A = secciones_piso[j-1][1]
        elements.append(Elements(nodes[i+(j-1)*nodos_ancho], nodes[i+j*nodos_ancho], A=A))

    if losas[j-1] == 1:
        #Defino elementos horizontales que conectan los ultimos nodos creados
        for i in range(nodos_ancho-1):
            nodos_actuales = len(nodes)
            elements.append(Elements(nodes[nodos_actuales-(nodos_ancho) + i], nodes[nodos_actuales-nodos_ancho + i + 1 ], q=-1000))

    base_v = vertical
#Ahora tomo la masa total de la estructra

M = Assembly(nodes, elements)

Des = Desplazamientos(nodes, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

import matplotlib.cm as cm
import matplotlib.colors as mcolors

def plot_original_structure_all_forces(nodes, elements):
    fig, ax = plt.subplots(figsize=(8, 6))

    # Obtener todas las áreas para normalización de colores
    areas = np.array([element.A for element in elements])
    norm = mcolors.Normalize(vmin=np.min(areas), vmax=np.max(areas))
    cmap = cm.viridis  # Puedes cambiar a otro mapa de colores como 'plasma', 'coolwarm', etc.

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
        if fy != 0:
            ax.quiver(x, y, 0, -1, angles='xy', scale_units='xy', scale=1, color='r',
                      label="Fuerza en Y" if node == nodes[0] else "")

        # Graficar momento como una flecha curva
        if m != 0:
            arc = plt.Circle((x, y), 0.2, color='purple', fill=False, linestyle='dashed')
            ax.add_patch(arc)
            ax.text(x + 0.25, y + 0.25, f'M={m:.2f}', fontsize=10, color='purple')

    # Barra de color para indicar valores de área
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    #cbar = plt.colorbar(sm, ax=ax)
    #cbar.set_label("Área de la sección [mm²]")

    # Configuración del gráfico
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title("Estructura Original con Colores según Área de la Sección")
    #ax.legend()
    ax.axis("equal")
    plt.grid(True)
    plt.show()

# Ejecutar la función para graficar la estructura con todas las fuerzas aplicadas
plot_original_structure_all_forces(nodes, elements)

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

# Ejecutar la función para graficar la estructura desplazada con escala ajustable
plot_deformed_structure(nodes, elements, scale=1000000)






