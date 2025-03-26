import numpy as np
from solver import Solver
from matrix_assembly import Assembly
from structure_assembly import Structure
from graph import plot_original_structure_all_forces
import matplotlib.pyplot as plt

structure_1 = Structure("Estructura 1")

#plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)

M = Assembly(structure_1.nodes, structure_1.elements)

Des = Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)



#Primero extraigo las deformaciones por nodos

contador = 0  # índice para self.uf_v

uf_v = Des.uf_v

for node in structure_1.nodes:
    for i, b in enumerate(node.boundary):
        if b == 0:  # Solo si el GDL está libre
            if i == 0 or i == 1:
                node.def_vector[i] += uf_v[contador] /1000
            else:
                node.def_vector[i] += uf_v[contador]  
            contador += 1

#Ahora genero los distintos vecotres deformacion

for element in structure_1.elements:
    #Con eso se cuales nodos conectan el elemento
    #Por lo tanto, puedo ensamblar el vector del elemento

    u = np.zeros(6)

    u[0:3] = element.n1.def_vector
    u[3:6] = element.n2.def_vector

    #Por lo que ahora genero los vecotres desplazamiento

    element.extractDisplacements(u)

    #element.plotGeometry(text=False, nodes=True, nodes_labels=False)



    #element.plotGeometry(text=False, nodes=True, nodes_labels=False)


import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

def plotStructure(beams, text=False, nodes=True, nodes_labels=False, deformada=True, escala=100000):
    """
    Función para graficar toda la estructura con todos los elementos (beams),
    mostrando la estructura original en un subgráfico y solo la deformada en otro subgráfico.
    
    :param beams: Lista de objetos `Frame2DOffsets` que representan los beams de la estructura.
    :param text: Indica si se debe mostrar el texto con la longitud útil.
    :param nodes: Indica si se deben mostrar los nodos en el gráfico.
    :param nodes_labels: Indica si se deben mostrar las etiquetas de los nodos.
    :param deformada: Indica si se debe mostrar la deformación de los elementos.
    :param escala: Factor de escala para la deformación.
    """
    if not beams:
        raise ValueError("La lista de beams no puede estar vacía.")
    
    # Crear los subgráficos
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))  # Dos subgráficos: uno para la estructura original, otro para la deformada
    
    # Inicializar las leyendas
    ax1_legend = []
    ax2_legend = []

    # Estilos comunes para los gráficos
    offset_color = 'orange'  # Color para los offsets rígidos
    elemento_color = 'blue'  # Color para el elemento útil
    deformada_color = 'r'  # Color para la deformada
    linewidth = 1  # Ancho de las líneas

    # Recorrer cada beam y graficarlo
    for beam in beams:
        # Coordenadas reales (nodos)
        xi_real, yi_real = beam.coord_i
        xf_real, yf_real = beam.coord_f

        # Coordenadas con offset (barra útil)
        xi, yi = beam.coord_i_offset
        xf, yf = beam.coord_f_offset
        
        # 2. Offsets rígidos (para la estructura original)
        ax1.plot([xi_real, xi], [yi_real, yi], offset_color, linewidth=linewidth)
        ax1.plot([xf, xf_real], [yf, yf_real], offset_color, linewidth=linewidth)

        # 3. Elemento útil (sin deformación)
        ax1.plot([xi, xf], [yi, yf], elemento_color , linewidth=linewidth)

        # 4. Nodos reales
        

        # 5. Etiquetas de nodo
        if nodes_labels:
            ax1.text(xi_real, yi_real, f"i ({xi_real:.2f}, {yi_real:.2f})", fontsize=9, ha='right', va='bottom')
            ax1.text(xf_real, yf_real, f"f ({xf_real:.2f}, {yf_real:.2f})", fontsize=9, ha='left', va='top')

        # 6. Texto del largo
        if text:
            xm = (xi + xf) / 2
            ym = (yi + yf) / 2
            ax1.text(xm, ym, f"L útil = {beam.L:.2f} m", fontsize=10, color=elemento_color, ha='center')

        # 7. Dibujar deformada interpolada
        if deformada and beam.u_global is not None:
            u_global = beam.u_global.flatten()

            # Transformar a coordenadas locales
            u_local = beam.tlg @ (beam.ts @ u_global)

            # Extraer grados de libertad locales
            ui_x, ui_y, theta_i, uj_x, uj_y, theta_j = u_local

            # Coordenadas locales a lo largo del eje
            L = beam.L
            n_points = 50
            x_vals = np.linspace(0, L, n_points)

            # Interpolación en X (axial)
            u_x_vals = (1 - x_vals / L) * ui_x + (x_vals / L) * uj_x

            # Funciones de forma para flexión (en Y)
            def N1(x): return 1 - 3*(x/L)**2 + 2*(x/L)**3
            def N2(x): return x*(1 - x/L)**2
            def N3(x): return 3*(x/L)**2 - 2*(x/L)**3
            def N4(x): return -x*(x/L)*(1 - x/L)

            v_vals = [N1(x)*ui_y + N2(x)*theta_i*L + N3(x)*uj_y + N4(x)*theta_j*L for x in x_vals]

            # Puntos locales deformados (X + Y interpolado)
            points_local = np.vstack([u_x_vals, v_vals]) * escala + np.vstack([x_vals, np.zeros_like(x_vals)])

            # Transformar a coordenadas globales
            points_global = beam.R @ points_local

            # Trasladar al sistema global real
            x0, y0 = beam.coord_i_offset
            x_def = points_global[0, :] + x0
            y_def = points_global[1, :] + y0

            ax2.plot(x_def, y_def, deformada_color + '-', linewidth=linewidth)

            # Nodo i
            p0i, p1i = beam.offset_rigido_deformado(beam.coord_i, u_global[0:3], beam.offset_i_global, escala)
            ax2.plot([p0i[0], p1i[0]], [p0i[1], p1i[1]], offset_color, linewidth=linewidth)

            # Nodo j
            p0j, p1j = beam.offset_rigido_deformado(beam.coord_f, u_global[3:6], beam.offset_j_global, escala)
            ax2.plot([p0j[0], p1j[0]], [p0j[1], p1j[1]], offset_color, linewidth=linewidth)

        #Dibujar los cachos rigidos deformados


    # Ajustes para gráficos
    ax1.set_aspect('equal')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.grid(True)
    
    ax2.set_aspect('equal')
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.grid(True)

    # Leyenda unificada para toda la estructura
    ax1.legend(['Offset rígido', 'Elemento útil'], loc='upper left', bbox_to_anchor=(1, 1))
    ax2.legend(['Deformada (x{})'.format(escala), 'Offset rígido'], loc='upper left', bbox_to_anchor=(1, 1))

    plt.tight_layout()
    plt.show()



plotStructure(structure_1.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=1000)