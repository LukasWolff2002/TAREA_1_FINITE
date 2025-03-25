import numpy as np
import matplotlib.pyplot as plt

fc = 28
E = 4700 * np.sqrt(fc)

gamma = 7800 #Kg/m3

class Elements:

    def __init__ (self, n1, n2, A=[200, 200], q=0, dxdy = [0.0, 0.0], E=E):
        #De base defino un area muy grande para las secciones que son axialmente rigidas

        self.n1 = n1
        self.n2 = n2
        self.coord_i = n1.coord
        self.coord_f = n2.coord
        self.q = q
        self.E = E
        self.A = A[0]*A[1]
        self.I = A[0]*A[1]**3/12
        self.dx = dxdy[0]
        self.dy = dxdy[1]
        self.L, self.angle = self.geometry()

        self.tbl = self.basicLocalTransformation()
        self.tlg = self.localGlobalTransformation()

        self.kb = self.basic_matrix()
        self.kl = self.local_matrix()
        self.k_global = self.global_matrix()
        self.ts = self.transformationStiffnessMatrix()
        self.ks_global = self.global_ts_matrix()
        self.Estructure_1()
        
    

    def geometry (self):
        coord_i_real = self.coord_i.copy()
        coord_f_real = self.coord_f.copy()

        delta = coord_f_real - coord_i_real
        length_total = np.linalg.norm(delta)
        direction_unit = delta / length_total
        angle = np.arctan2(delta[1], delta[0])

        # Matriz de transformación local → global
        c = direction_unit[0]
        s = direction_unit[1]
        R = np.array([[c, -s],
                      [s,  c]])

        # Offsets locales independientes para cada nodo
        offset_i_local = np.array([self.dx, self.dy])      # nodo i
        offset_j_local = np.array([-self.dx, self.dy])     # nodo j (inverso en X local)

        # Transformarlos a coordenadas globales
        offset_i_global = R @ offset_i_local
        offset_j_global = R @ offset_j_local

        self.offset_i_global = offset_i_global
        self.offset_j_global = offset_j_global

        self.L_offset_i = np.linalg.norm(offset_i_global)
        self.L_offset_j = np.linalg.norm(offset_j_global)



        # Aplicar a los extremos
        coord_i_offset = coord_i_real + offset_i_global
        coord_f_offset = coord_f_real + offset_j_global

        # Guardar valores
        self.coord_i_offset = coord_i_offset
        self.coord_f_offset = coord_f_offset
        self.R = R
        self.angle = angle

        length_effective = np.linalg.norm(coord_f_offset - coord_i_offset)
        return length_effective, angle

    def basic_matrix(self):
        L = self.L * 1000
        A = self.A
        E = self.E
        I = self.I
        
        Kb = np.array([[A*E/L, 0, 0],
                       [0, 4*E*I/L, 2*E*I/L],
                       [0, 2*E*I/L, 4*E*I/L]])
        
        return Kb
    
    def basicLocalTransformation(self):
        L = self.L * 1000
        Tbl = np.array([
            [-1, 0, 0, 1, 0, 0],
            [0, 1/L, 1, 0, -1/L, 0],
            [0, 1/L, 0, 0, -1/L, 1]
            ]) 
        return Tbl
    
    def localGlobalTransformation(self):
        c = np.cos(self.angle)
        s = np.sin(self.angle)

        Tlg = np.array([[ c, s, 0,  0, 0, 0],
                        [-s, c, 0,  0, 0, 0],
                        [ 0, 0, 1,  0, 0, 0],
                        [ 0, 0, 0,  c, s, 0],
                        [ 0, 0, 0, -s, c, 0],
                        [ 0, 0, 0,  0, 0, 1]
                        ])
        return Tlg

    def local_matrix (self):
    
        Tbl = self.tbl
        Kb = self.kb
        Kl = Tbl.T @ Kb @ Tbl
        
        return Kl
    

    def transformationStiffnessMatrix(self):
        # Offsets en coordenadas locales
        offset_i_local = np.array([self.dx, self.dy])     # Nodo i
        offset_j_local = np.array([-self.dx, self.dy])    # Nodo j (inverso en X local)

        # Transformación local → global
        offset_i_global = self.R @ offset_i_local
        offset_j_global = self.R @ offset_j_local

        dxi, dyi = offset_i_global
        dxj, dyj = offset_j_global

        Ts = np.array([
            [1, 0, -dyi, 0, 0,   0],
            [0, 1,  dxi, 0, 0,   0],
            [0, 0,    1, 0, 0,   0],
            [0, 0,    0, 1, 0, -dyj],  # ← sin cambio de signo
            [0, 0,    0, 0, 1,  dxj],  # ← sin cambio de signo
            [0, 0,    0, 0, 0,    1]
        ])

        return Ts

        
    def global_matrix (self):
        k = self.kl
        tlg = self.tlg
        Kg = tlg @ k @ tlg.T
        
        return Kg
    
    def global_ts_matrix (self):
        kg = self.k_global
        ts = self.ts
        
        Ks = ts.T @ kg @ ts
        
        return Ks
    
                      
    
    def Estructure_1(self):
        if self.q != 0:

            m = (self.q*self.L**2)/12 #
            v = (self.q*self.L)/2 

            #Ahora debo encontrar cual es el nodo de la izquierda y el de la derecha
            if self.n1.coord[0] < self.n2.coord[0]:
                
                self.n1.force_vector = self.n1.force_vector + np.array([0, v, m])
                self.n2.force_vector = self.n2.force_vector + np.array([0, v, -m])
            
            else:
                self.n1.force_vector += np.array([0, v, -m])
                self.n2.force_vector += np.array([0, v, m])

    def offset_rigido_deformado(self, nodo_real, u_nodo, offset_global):

        ux, uy, theta = u_nodo
        p_deformado = nodo_real + np.array([ux, uy])

        # Rotación del offset alrededor del nodo
        R_theta = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta),  np.cos(theta)]
        ])
        offset_rotado = R_theta @ offset_global

        p_final = p_deformado + offset_rotado
        return p_deformado, p_final


    def extractDisplacements (self, u):
        #Aqui el vector u es el vector de desplazamientos de los nodos
      
        u_global = u
        
        Tlg = self.tlg
        u_local = Tlg @ u_global
        Tbl = self.tbl
        u_basic = Tbl @ u_local

        self.u_global = u_global
        self.u_local = u_local
        self.u_basic = u_basic

        # Calcular los extremos deformados del elemento útil
        _, p1i = self.offset_rigido_deformado(
            self.coord_i, u[0:3], self.offset_i_global
        )
        _, p1j = self.offset_rigido_deformado(
            self.coord_f, u[3:6], self.offset_j_global
        )

        # Guardar vector deformado como array de 2x2: [[xi, yi], [xf, yf]]
        self.u_corrected = np.array([p1i, p1j])

    def plotGeometry(self, ax=None, text=False, nodes=True, nodes_labels=False, deformada=True, escala=100000):

        if ax is None:
            fig, ax = plt.subplots()

        # Coordenadas reales (nodos)
        xi_real, yi_real = self.coord_i
        xf_real, yf_real = self.coord_f

        # Coordenadas con offset (barra útil)
        xi, yi = self.coord_i_offset
        xf, yf = self.coord_f_offset

        # 1. Eje total entre nodos reales
        ax.plot([xi_real, xf_real], [yi_real, yf_real], 'k--', linewidth=1)

        # 2. Offsets rígidos
        ax.plot([xi_real, xi], [yi_real, yi], color='orange', linewidth=3)
        ax.plot([xf, xf_real], [yf, yf_real], color='orange', linewidth=3)

        # 3. Elemento útil (sin deformación)
        ax.plot([xi, xf], [yi, yf], 'b-', linewidth=2)

        # 4. Nodos reales
        if nodes:
            ax.plot(xi_real, yi_real, 'ro')
            ax.plot(xf_real, yf_real, 'ro')

        # 5. Etiquetas de nodo
        if nodes_labels:
            ax.text(xi_real, yi_real, f"i ({xi_real:.2f}, {yi_real:.2f})", fontsize=9, ha='right', va='bottom')
            ax.text(xf_real, yf_real, f"f ({xf_real:.2f}, {yf_real:.2f})", fontsize=9, ha='left', va='top')

        # 6. Texto del largo
        if text:
            xm = (xi + xf) / 2
            ym = (yi + yf) / 2
            ax.text(xm, ym, f"L útil = {self.L:.2f} m", fontsize=10, color='darkgreen', ha='center')

            # 7. Dibujar deformada interpolada
        if deformada and self.u_global is not None:
            u_global = self.u_global.flatten()

            # Transformar a coordenadas locales
            u_local = self.tlg @ (self.ts @ u_global)

            # Extraer grados de libertad locales
            ui_x, ui_y, theta_i, uj_x, uj_y, theta_j = u_local

            # Coordenadas locales a lo largo del eje
            L = self.L
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
            points_global = self.R @ points_local

            # Trasladar al sistema global real
            x0, y0 = self.coord_i_offset
            x_def = points_global[0, :] + x0
            y_def = points_global[1, :] + y0

            ax.plot(x_def, y_def, 'r--', linewidth=2)

            # Escala
            escala = escala

            # Nodo i
            p0i, p1i = self.offset_rigido_deformado(self.coord_i, u_global[0:3]*escala, self.offset_i_global*escala)
            ax.plot([p0i[0], p1i[0]], [p0i[1], p1i[1]], 'r--', linewidth=3)

            # Nodo j
            p0j, p1j = self.offset_rigido_deformado(self.coord_f, u_global[3:6]*escala, self.offset_j_global*escala)
            ax.plot([p0j[0], p1j[0]], [p0j[1], p1j[1]], 'r--', linewidth=3)

        ax.set_aspect('equal')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.grid(True)
        # Mueve la leyenda fuera del gráfico, al lado derecho
        #ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
 
        
        #plt.show()



                


    
    
    
    


    