
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

def plot_deformed_structures(structure_list, text=False, nodes=True, nodes_labels=False, deformada=True, Escala=[1,1,1]):
    
    beams_list = []

    for structure in structure_list:
        beams_list.append(structure.elements)
        
    # Crear los subgráficos
    fig, axes = plt.subplots(1, len(beams_list), figsize=(18, 6))  # Ajustar tamaño basado en número de estructuras

    # Si solo hay una estructura, `axes` es un solo objeto, así que lo convertimos en una lista
    if len(beams_list) == 1:
        axes = [axes]

    # Estilos comunes para los gráficos
    offset_color = 'orange'  # Color para los offsets rígidos
    elemento_color = 'blue'  # Color para el elemento útil
    deformada_color = 'r'  # Color para la deformada
    linewidth = 1  # Ancho de las líneas

    for i, beams in enumerate(beams_list):

        ax = axes[i]  # Seleccionar el gráfico correspondiente para esta estructura

        escala = Escala[i]
        # Recorrer cada beam y graficarlo
        for beam in beams:
            # Coordenadas reales (nodos)
            xi_real, yi_real = beam.coord_i
            xf_real, yf_real = beam.coord_f

            # Coordenadas con offset (barra útil)
            xi, yi = beam.coord_i_offset
            xf, yf = beam.coord_f_offset
            
            # 7. Dibujar deformada interpolada
            if deformada and beam.u_global is not None:
                u_global = beam.u_global.flatten()

                # Transformar a coordenadas locales
                u_local = beam.tlg @ (beam.ts @ u_global)

                # Extraer grados de libertad locales
                ui_x, ui_y, theta_i, uj_x, uj_y, theta_j = u_local

                # Coordenadas locales a lo largo del eje
                L = beam.L

                theta_i = theta_i / L
                theta_j = theta_j / L

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

                ax.plot(x_def, y_def, deformada_color + '-', linewidth=linewidth)
                u_nodo_scaled = u_global[0:3]
                u_nodo_scaled[:3] *= escala  # solo desplazamientos, no rotación

                # Nodo i
                p0i, p1i = beam.offset_rigido_deformado(beam.coord_i, u_nodo_scaled, beam.offset_i_global, escala)
                ax.plot([p0i[0], p1i[0]], [p0i[1], p1i[1]], offset_color, linewidth=linewidth)

                u_nodo_scaled = u_global[3:]
                u_nodo_scaled[:3] *= escala  # solo desplazamientos, no rotación

                # Nodo j
                p0j, p1j = beam.offset_rigido_deformado(beam.coord_f, u_nodo_scaled, beam.offset_j_global, escala)
                ax.plot([p0j[0], p1j[0]], [p0j[1], p1j[1]], offset_color, linewidth=linewidth)

            #Dibujar los cachos rigidos deformados


        # Ajustes para gráficos
        ax.set_aspect('equal')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        titulo = structure_list[i].nombre
        ax.set_title(f"{str(titulo)} - Escala: {escala}:1")
        ax.grid(False)

        #ax.legend(['Deformada (x{})'.format(escala), 'Offset rígido'], loc='upper left', bbox_to_anchor=(1, 1))
    plt.suptitle("Estructuras deformadas",fontsize=20)
    plt.tight_layout()
    plt.show()


def plotStructureWithMomentDiagram(structure_list, Escala, show_structure):
    
    fig, axes = plt.subplots(1, len(structure_list), figsize=(6 * len(structure_list), 6))

    if len(structure_list) == 1:
        axes = [axes]

    for i, structure in enumerate(structure_list):
        ax = axes[i]
        escala_m = Escala[i] if i < len(Escala) else Escala[-1]

        for beam in structure.elements:
            xi, yi = beam.coord_i_offset
            xf, yf = beam.coord_f_offset
            L = beam.L

            if not hasattr(beam, "f_local"):
                beam.extractForces()

            f_local = beam.f_local

            M1, M2 = f_local[2], f_local[5]

            # Determinar si es horizontal o vertical
            is_horizontal = np.isclose(yi, yf, atol=1e-6)

            x_local = np.linspace(0, L, 100)

            if is_horizontal and beam.q != 0:
                q = beam.q
                # Momentos parabólicos para vigas horizontales con carga distribuida
                M_vals = M1 * (1 - x_local/L) + M2 * (x_local/L) - (q * x_local * (L - x_local)) / 2
            else:
                # Momentos lineales para columnas y vigas sin carga distribuida
                M_vals = np.linspace(M1, M2, len(x_local))

            # Escalar el momento
            y_local = -M_vals * escala_m
            points_local = np.vstack([x_local, y_local])
            points_global = beam.R @ points_local
            x_global = points_global[0, :] + xi
            y_global = points_global[1, :] + yi

            # Graficar momento flectores
            ax.plot(x_global, y_global, 'purple', linewidth=1, label='Momento' if beam == structure.elements[0] else "")
            
            # Dibujar estructura original si se requiere
            if show_structure:
                ax.plot([xi, xf], [yi, yf], 'b-', linewidth=0.5)

        ax.set_aspect('equal')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        titulo = getattr(structure, "nombre", f"Estructura {i+1}")
        ax.set_title(f"{titulo} - Escala: {escala_m}:1")

        ax.grid(False)
        ax.legend()

    plt.suptitle("Diagramas de momentos", fontsize=20)
    plt.tight_layout()
    plt.show()



def plotStructureWithShearDiagram(structure_list, Escala=[1e-4,1e-4,1e-4], show_structure=True):
    
    beams_list = []

    for structure in structure_list:
        beams_list.append(structure.elements)
        
    # Crear los subgráficos
    fig, axes = plt.subplots(1, len(beams_list), figsize=(18, 6))  # Ajustar tamaño basado en número de estructuras

    # Si solo hay una estructura, `axes` es un solo objeto, así que lo convertimos en una lista
    if len(beams_list) == 1:
        axes = [axes]

    for i, beams in enumerate(beams_list):

        ax = axes[i]  # Seleccionar el gráfico correspondiente para esta estructura

        escala_corte = Escala[i]
        for beam in beams:
            # Extraer extremos de la barra útil
            xi, yi = beam.coord_i_offset
            xf, yf = beam.coord_f_offset
            L = beam.L

            # Asegurar que tiene fuerzas internas
            if not hasattr(beam, "f_local"):
                beam.extractForces()

            if xi == xf:
                #Tengo un elemento vertical

                # Fuerzas cortantes locales
                V1 = beam.f_local[1] 
                V2 = -beam.f_local[4] 
            
            else:
                # Fuerzas cortantes locales
                V1 = beam.f_local[0]
                V2 = beam.f_local[3]

            # Coordenadas locales a lo largo del eje
            x_local = np.linspace(0, L, 50)
            V_vals = (1 - x_local / L) * V1 + (x_local / L) * V2  # interpolación lineal

            # Coordenadas locales (x, y corte)
            points_local = np.vstack([x_local, -V_vals * escala_corte])  # signo negativo: convención hacia arriba

            # Transformar a global
            points_global = beam.R @ points_local
            x_global = points_global[0, :] + xi
            y_global = points_global[1, :] + yi

            # Dibujar la línea del diagrama
            ax.plot(x_global, y_global, color='green', linewidth=1, label='Corte' if beam == beams[0] else "")

            # Dibujar la estructura si se desea
            if show_structure:
                ax.plot([xi, xf], [yi, yf], 'b-', linewidth=0.5)
                

            # Marcar valores de corte al inicio y final
            #ax.text(x_global[0], y_global[0], f'{V1:.1f}', fontsize=8, color='darkgreen', ha='right')
            #ax.text(x_global[-1], y_global[-1], f'{V2:.1f}', fontsize=8, color='darkgreen', ha='left')

        ax.set_aspect('equal')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f"{structure_list[i].nombre} - Escala: {escala_corte}:1")
        ax.grid(False)
        ax.legend()
    
    plt.suptitle("Diagramas de corte",fontsize=20)
    plt.tight_layout()
    plt.show()

def plotStructureWithAxialDiagram(structure_list, Escala=[1,1,1], show_structure=True):
    
    beams_list = []

    for structure in structure_list:
        beams_list.append(structure.elements)
        
    # Crear los subgráficos
    fig, axes = plt.subplots(1, len(beams_list), figsize=(18, 6))  # Ajustar tamaño basado en número de estructuras

    # Si solo hay una estructura, `axes` es un solo objeto, así que lo convertimos en una lista
    if len(beams_list) == 1:
        axes = [axes]

    for i, beams in enumerate(beams_list):

        ax = axes[i]  # Seleccionar el gráfico correspondiente para esta estructura

        escala_axial = Escala[i]

        for beam in beams:
            # Extremos útiles
            xi, yi = beam.coord_i_offset
            xf, yf = beam.coord_f_offset
            L = beam.L

            # Asegurarse de tener fuerzas internas
            if not hasattr(beam, "f_basic"):
                beam.extractForces()

            N = beam.f_basic[0]  # esfuerzo axial

            # Coordenadas locales del rectángulo de esfuerzo axial
            x_local = np.array([0, L])
            y_local = np.array([-N * escala_axial, -N * escala_axial])

            # Crear puntos del diagrama
            points_local = np.vstack([x_local, y_local])
            points_global = beam.R @ points_local
            x_global = points_global[0, :] + xi
            y_global = points_global[1, :] + yi

            # Dibujar línea de esfuerzo axial
            ax.plot(x_global, y_global, color='brown', linewidth=1, label='Axial' if beam == beams[0] else "")

            if show_structure:
                ax.plot([xi, xf], [yi, yf], 'b-', linewidth=0.5)
                #ax.plot(*beam.coord_i, 'ro')
                #ax.plot(*beam.coord_f, 'ro')

            # Mostrar valor
            #ax.text(x_global[0], y_global[0], f'{N:.1f}', fontsize=8, color='saddlebrown', ha='right')
            #ax.text(x_global[1], y_global[1], f'{N:.1f}', fontsize=8, color='saddlebrown', ha='left')
        
        ax.set_aspect('equal')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.grid(False)
        ax.set_title(f"{structure_list[i].nombre} - Escala: {escala_axial}:1")
        ax.legend()

    plt.suptitle("Diagramas de esfuerzos axiales",fontsize=20)
    plt.tight_layout()
    plt.show()

