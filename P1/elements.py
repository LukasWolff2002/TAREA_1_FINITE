import numpy as np

fc = 28
E = 4700 * np.sqrt(fc)

gamma = 7800 #Kg/m3

class Elements:

    def __init__ (self, n1, n2, A=[100000, 100000], q=0):
        #De base defino un area muy grande para las secciones que son axialmente rigidas

        self.n1 = n1
        self.n2 = n2
        self.coords_i = n1.coord
        self.coords_f = n2.coord
        self.q = q
        #self.E = E
        self.A = A[0]*A[1]
        self.I = A[0]*A[1]**3/12
        #self.L = L

        self.L, self.angle = self.length_angle()
        self.k_local = self.local_matrix()
        self.k_global = self.global_matrix()
        self.Estructure_1()

    def length_angle (self):
        length = (np.linalg.norm(self.coords_f - self.coords_i)) 
        angle = np.arctan2(self.coords_f[1] - self.coords_i[1], self.coords_f[0] - self.coords_i[0])
        return length, angle

    def local_matrix (self):
        L = self.L * 1000
        I = self.I
        A = self.A
        
        
        #Por lo tanto puedo definir mi matriz K en coordenadas globales
        k = np.array([[(A*E)/L, 0, 0, -(A*E)/L, 0, 0],
                      [0, (12*E*I)/(L**3), (6*E*I)/(L**2), 0, -(12*E*I)/(L**3), (6*E*I)/(L**2)],
                      [0, (6*E*I)/(L**2), (4*E*I)/L, 0, -(6*E*I)/(L**2), (2*E*I)/L],
                      [-(A*E)/L, 0, 0, (A*E)/L, 0, 0],
                      [0, -(12*E*I)/(L**3), -(6*E*I)/(L**2), 0, (12*E*I)/(L**3), -(6*E*I)/(L**2)],
                      [0, (6*E*I)/(L**2), (2*E*I)/L, 0, -(6*E*I)/(L**2), (4*E*I)/L]])
        
        return k
        
    def global_matrix (self):
        k = self.k_local

        c = np.cos(self.angle)
        s = np.sin(self.angle)

        T = np.array([[c, -s, 0, 0, 0, 0], 
                      [s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, -s, 0],
                      [0, 0, 0, s, c, 0],
                      [0, 0, 0, 0, 0, 1]])
        
        Ke = T @ k @ T.T
        
        return Ke
    
    def Estructure_1(self):
        if self.q != 0:

            m = (self.q*self.L**2)/12
            v = (self.q*self.L)/2

            #Ahora debo encontrar cual es el nodo de la izquierda y el de la derecha
            if self.n1.coord[0] < self.n2.coord[0]:
                
                self.n1.force_vector = self.n1.force_vector + np.array([0, v, m])
                self.n2.force_vector = self.n2.force_vector + np.array([0, v, -m])
            
            else:
                self.n1.force_vector += np.array([0, v, -m])
                self.n2.force_vector += np.array([0, v, m])

        
 