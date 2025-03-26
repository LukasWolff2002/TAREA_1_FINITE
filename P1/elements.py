import numpy as np
import matplotlib.pyplot as plt

fc = 28
E = 4700 * np.sqrt(fc)

gamma = 7800 #Kg/m3

class Elements:

    def __init__ (self, n1, n2, AI, q=0, dxdy = [0.0, 0.0], E=E):
        #De base defino un area muy grande para las secciones que son axialmente rigidas

        self.n1 = n1
        self.n2 = n2
        self.coord_i = n1.coord
        self.coord_f = n2.coord
        self.q = q
        self.E = E
        self.A = AI[0]
        self.I = AI[1]
        self.dx = dxdy[0]
        self.dy = dxdy[1]
        self.L, self.angle, self.Peso = self.geometry()

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
        Peso = -length_total * (self.A/1000000) * gamma * 9.81
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
        return length_effective, angle, Peso

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

        if self.n1.coord[0] == self.n2.coord[0]:

            #Es una columna
            #Es una viga
            q = self.Peso

            n = ((q*self.L)/2)
       
            self.n1.force_vector = self.n1.force_vector + np.array([0, n, 0])
            self.n2.force_vector = self.n2.force_vector + np.array([0, n, 0])
            
        
        else:

           #Es una viga
           q = self.q + self.Peso
           m = (q*self.L**2)/12 #
           v = (q*self.L)/2 
           #Ahora debo encontrar cual es el nodo de la izquierda y el de la derecha
           if self.n1.coord[0] < self.n2.coord[0]:
               
               self.n1.force_vector = self.n1.force_vector + np.array([0, v, m])
               self.n2.force_vector = self.n2.force_vector + np.array([0, v, -m])
           
           else:
               self.n1.force_vector += np.array([0, v, -m])
               self.n2.force_vector += np.array([0, v, m])

    def offset_rigido_deformado(self, nodo_real, u_nodo, offset_global, escala=1):

        ux, uy, theta = u_nodo
        p_deformado = nodo_real + np.array([ux, uy]) 

        # Rotación del offset alrededor del nodo
        R_theta = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta),  np.cos(theta)]
        ])
        offset_rotado = (R_theta @ offset_global) 
     

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

   


                


    
    
    
    


    