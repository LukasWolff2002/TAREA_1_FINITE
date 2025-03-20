import numpy as np

fc = 28
E = 4700 * np.sqrt(fc)
A = 300 * 300
I = 300 * 300 ** 3 / 12

A = A * 1000

class Elements:

    def __init__ (self, n1, n2):
        self.n1 = n1.n
        self.n2 = n2.n
        self.coords_i = n1.coord
        self.coords_f = n2.coord
        #self.E = E
        #self.A = A
        #self.I = I
        #self.L = L

        self.L, self.angle = self.length_angle()
        self.k_local = self.local_matrix()
        self.k_global = self.global_matrix()

    def length_angle (self):
        length = (np.linalg.norm(self.coords_f - self.coords_i)) * 1000
        angle = np.arctan2(self.coords_f[1] - self.coords_i[1], self.coords_f[0] - self.coords_i[0])
        return length, angle

    def local_matrix (self):
        #L = self.length
        
        
        #Por lo tanto puedo definir mi matriz K en coordenadas globales
        k = np.array([[(A*E)/self.L, 0, 0, -(A*E)/self.L, 0, 0],
                      [0, (12*E*I)/(self.L**3), (6*E*I)/(self.L**2), 0, -(12*E*I)/(self.L**3), (6*E*I)/(self.L**2)],
                      [0, (6*E*I)/(self.L**2), (4*E*I)/self.L, 0, -(6*E*I)/(self.L**2), (2*E*I)/self.L],
                      [-(A*E)/self.L, 0, 0, (A*E)/self.L, 0, 0],
                      [0, -(12*E*I)/(self.L**3), -(6*E*I)/(self.L**2), 0, (12*E*I)/(self.L**3), -(6*E*I)/(self.L**2)],
                      [0, (6*E*I)/(self.L**2), (2*E*I)/self.L, 0, -(6*E*I)/(self.L**2), (4*E*I)/self.L]])
        
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
 