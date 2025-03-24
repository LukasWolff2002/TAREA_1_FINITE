import numpy as np

fc = 28
E = 4700 * np.sqrt(fc)

gamma = 7800 #Kg/m3

class Elements:

    def __init__ (self, n1, n2, A=[200, 200], q=0, dxdy = [0.0, 0.0], E=E):
        #De base defino un area muy grande para las secciones que son axialmente rigidas

        self.n1 = n1
        self.n2 = n2
        self.coords_i = n1.coord
        self.coords_f = n2.coord
        self.q = q
        self.E = E
        self.A = A[0]*A[1]
        self.I = A[0]*A[1]**3/12
        self.dx = dxdy[0]
        self.dy = dxdy[1]
        self.L, self.angle = self.length_angle()

        self.tbl = self.basicLocalTransformation()
        self.tlg = self.localGlobalTransformation()

        self.kb = self.basic_matrix()
        self.kl = self.local_matrix()
        self.tgo = self.tgo_matrix()
        self.k_global = self.global_matrix()
        self.k_global_tgo = self.global_tgo_matrix()
        self.Estructure_1()
        
        
        
        
        

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

        
    def global_matrix (self):
        k = self.kl
        tlg = self.tlg
        
        Kg = tlg @ k @ tlg.T
        
        return Kg
    
    def global_tgo_matrix (self):
        k = self.kl
        tgo = self.tgo
        tlg = self.tlg
        
        Ke = tlg @ k @ tlg.T
        Ke = tgo.T @ Ke @ tgo
        
        return Ke
    
                      
    
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

                


    
    
    
    


    