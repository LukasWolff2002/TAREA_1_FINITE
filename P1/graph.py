
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

def plot_deformed_structure(beams, text=False, nodes=True, nodes_labels=False, deformada=True, escala=100000):

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

            u_nodo_scaled = u_global[0:3]
            u_nodo_scaled[:2] *= escala  # solo desplazamientos, no rotación

            # Nodo i
            p0i, p1i = beam.offset_rigido_deformado(beam.coord_i, u_nodo_scaled, beam.offset_i_global, escala)
            ax2.plot([p0i[0], p1i[0]], [p0i[1], p1i[1]], offset_color, linewidth=linewidth)

            u_nodo_scaled = u_global[3:]
            u_nodo_scaled[:2] *= escala  # solo desplazamientos, no rotación

            # Nodo j
            p0j, p1j = beam.offset_rigido_deformado(beam.coord_f, u_nodo_scaled, beam.offset_j_global, escala)
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
