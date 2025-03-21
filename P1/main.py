import numpy as np
from nodes import Node
from elements import Elements
from desplazamientos import Desplazamientos
from assembly import Assembly
from graph import plot_original_structure_all_forces, plot_deformed_structure, plot_structure_with_local_displacements
import matplotlib.pyplot as plt

nodes = []
elements = []

espaciado_h = 9.14 #m
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

#Seccion (mm, mm)   I y D       Centrales   
secciones_piso = [[[300, 300], [100, 200]], #Piso 1
                  [[300, 300], [100, 200]], #Piso 2
                  [[300, 300], [100, 200]], #Piso 3
                  [[300, 300], [100, 200]], #Piso 4
                  [[300, 300], [100, 200]], #Piso 5
                  [[300, 300], [100, 200]], #Piso 6
                  [[300, 300], [100, 200]], #Piso 7
                  [[300, 300], [100, 200]], #Piso 8
                  [[300, 300], [100, 200]], #Piso 9
                  [[300, 300], [100, 200]], #Piso 10
                  [[300, 300], [100, 200]], #Piso 11
                  [[300, 300], [100, 200]], #Piso 12
                  [[300, 300], [100, 200]], #Piso 13
                  [[300, 300], [100, 200]]  #Piso 14
                  ]

cargas_q = [26,
            26,
            0,
            26,
            26,
            0,
            26,
            26,
            0,
            26,
            26,
            0,
            26,
            23.1]
            
            

nodos_ancho = 6
pisos = 14

base_v = 0
#Defino los nodos base
for i in range(nodos_ancho):
    nodes.append(Node(i, np.array([i*espaciado_h, 0]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))


#Defino los nodos de los pisos
for j in range(1, pisos+1):
    espaciado_v = Espaciado_v[j-1]
    vertical = base_v + espaciado_v

    if j == 1:
        for i in range(nodos_ancho):
            if i == 0:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

            elif i == nodos_ancho-1:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

            else:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

    else:
        for i in range(nodos_ancho):

            if i == 0:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

            else:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

    #Ahora conecto verticalmente los pisos
    for i in range(nodos_ancho):
        if i == 0:
            A = secciones_piso[j-1][0]
        elif i == nodos_ancho-1:
            A = secciones_piso[j-1][0]
        else:
            A = secciones_piso[j-1][1]
        elements.append(Elements(nodes[i+(j-1)*nodos_ancho], nodes[i+j*nodos_ancho], A=A))

    if losas[j-1] == 1:
        #Defino elementos horizontales que conectan los ultimos nodos creados
        for i in range(nodos_ancho-1):
            nodos_actuales = len(nodes)
            elements.append(Elements(nodes[nodos_actuales-(nodos_ancho) + i], nodes[nodos_actuales-nodos_ancho + i + 1 ], q=-1000*cargas_q[j-1]))

    base_v = vertical

#Ahora tomo la masa total de la estructra

M = Assembly(nodes, elements)

Des = Desplazamientos(nodes, elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

# Ejecutar la función para graficar la estructura c on todas las fuerzas aplicadas
#plot_original_structure_all_forces(nodes, elements)

# Ejecutar la función para graficar la estructura desplazada con escala ajustable
#plot_deformed_structure(nodes, elements, scale=50)

M = Assembly(nodes, elements, fixed = False)

Des = Desplazamientos(nodes, elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

#plot_deformed_structure(nodes, elements, scale=50)


# Ejecutar la función para graficar la estructura con las funciones de forma



#plot_structure_with_local_displacements(nodes, elements, scale=50)
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

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
        x_local = np.linspace(0, L, 100)  # Puntos a lo largo de la longitud del elemento

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



plot_deformed_structure(nodes, elements, scale=100)