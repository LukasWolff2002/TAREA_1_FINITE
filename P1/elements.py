import numpy as np

fc = 28
E = 4700 * np.sqrt(fc)

gamma = 7800 #Kg/m3

class Elements:

    def __init__ (self, n1, n2, A=[200, 200], q=0):
        #De base defino un area muy grande para las secciones que son axialmente rigidas

        self.n1 = n1
        self.n2 = n2
        self.coords_i = n1.coord
        self.coords_f = n2.coord
        self.q = q
        self.E = E
        self.A = A[0]*A[1]
        self.I = A[0]*A[1]**3/12
        #self.L = L

        self.L, self.angle = self.length_angle()
        self.kb = self.basic_matrix()
        self.k_local = self.local_matrix()
        self.tlg = self.localGlobalTransformation()
        self.tgo = self.tgo_matrix()
        self.k_global = self.global_matrix()
        self.k_global_tgo = self.global_tgo_matrix()
        self.Estructure_1()
        self.ug = np.array([0, 0, 0, 0, 0, 0])
        
        

    def length_angle (self):
        length = (np.linalg.norm(self.coords_f - self.coords_i)) 
        angle = np.arctan2(self.coords_f[1] - self.coords_i[1], self.coords_f[0] - self.coords_i[0])

        return length, angle

    def basic_matrix(self):
        L = self.L * 1000
        A = self.A
        E = self.E
        I = self.I
        
        Kb = np.array([[A*E/L, 0, 0],
                       [0, 4*E*I/L, 2*E*I/L],
                       [0, 2*E*I/L, 4*E*I/L]])
        
        return Kb

    def local_matrix (self):
        L = self.L * 1000
        I = self.I
        A = self.A
        
        Kb = self.kb
        
        Tbl = np.array([
            [-1, 0, 0, 1, 0, 0],
            [0, 1/L, 1, 0, -1/L, 0],
            [0, 1/L, 0, 0, -1/L, 1]
            ]) 
        
        Kl = Tbl.T @ Kb @ Tbl
        
        return Kl

    def tgo_matrix (self):
        dx = 1
        dy = 0

        tgo = np.array([[1, 0, -dy, 0, 0, 0],
                        [0, 1, dx, 0, 0, 0],
                        [0, 0, 1, 0, 0, 0],
                        [0, 0, 0, 1, 0, -dy],
                        [0, 0, 0, 0, 1, dx],
                        [0, 0, 0, 0, 0, 1]])
        
        return tgo
    
    def localGlobalTransformation(self):
        c = np.cos(self.angle)
        s = np.sin(self.angle)

        Tlg = np.array([
        [ c, s, 0,  0, 0, 0],
        [-s, c, 0,  0, 0, 0],
        [ 0, 0, 1,  0, 0, 0],
        [ 0, 0, 0,  c, s, 0],
        [ 0, 0, 0, -s, c, 0],
        [ 0, 0, 0,  0, 0, 1]
        ])
        return Tlg

        
    def global_matrix (self):
        k = self.k_local
        

        c = np.cos(self.angle)
        s = np.sin(self.angle)

        T = np.array([[c, s, 0, 0, 0, 0], 
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])
        
        Kg = T @ k @ T.T
        
        return Kg, T
    
    def global_tgo_matrix (self):
        k = self.k_local
        tgo = self.tgo

        c = np.cos(self.angle)
        s = np.sin(self.angle)

        T = np.array([[c, s, 0, 0, 0, 0], 
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])
        
        Ke = T @ k @ T.T
        Ke = tgo.T @ Ke @ tgo
        
        return Ke, T
    
                      
    
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


    def _extractDisplacements(self, u): #Extraer los desplazamientos de los nodos
        u_global = u[self.dof_indices].reshape((6, 1))  # 6 GDL por elemento
        Tlg = self.localGlobalTransformation()
        u_local = Tlg @ u_global
        Tbl = self.basicLocalTransformation()
        u_basic = Tbl @ u_local
        return u_global, u_local, u_basic
    
    def _calculateBasicForces(self, u):
        _, _, u_basic = self._extractDisplacements(u)
        Kb = self.stiffness_matrix_basic()
        f_basic = Kb @ u_basic
        return f_basic


    def _calculateLocalForces(self, u):
        _, u_local, _ = self._extractDisplacements(u)
        Kl = self.localStiffnessMatrix()
        f_local = Kl @ u_local
        return f_local

    def forceRecovery(self, u, printSummary=True):
        f_basic = self._calculateBasicForces(u)
        f_local = self._calculateLocalForces(u)
        Tlg = self.localGlobalTransformation()
        f_global = Tlg.T @ f_local

        if printSummary:
            print("Fuerzas básicas:\n", f_basic)
            print("Fuerzas locales:\n", f_local)
            print("Fuerzas globales:\n", f_global)

        return f_basic, f_local, f_global


    def local_displassments (self):
        n1 = self.n1
        n2 = self.n2

        c = np.cos(self.angle)
        s = np.sin(self.angle)

        ul = np.array([n1.def_vector[0], n1.def_vector[1], n1.def_vector[2], 
                       n2.def_vector[0], n2.def_vector[1], n2.def_vector[2]])

        T = np.array([[c, s, 0, 0, 0, 0], 
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])

        ug = T.T @ ul

        return ug

        
 