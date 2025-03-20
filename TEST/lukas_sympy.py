import numpy as np
import sympy as sp


E, A, I, L = sp.symbols('E A, I, L')

class Node:

    def __init__ (self, n, coord, ndof=3):
        self.n = n
        self.coord = coord
        self.ndof = ndof
        self.idx = self.index()

    def index (self):
        #Debo generar un array desde n hasta n+ndof, con np.array
        return np.linspace((self.n*self.ndof), (self.n*self.ndof)+(self.ndof-1), 3, dtype=int)
    
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

        self.length, self.angle = self.length_angle()
        self.k_local = self.local_matrix()
        self.k_global = self.global_matrix()

    def length_angle (self):
        length = np.linalg.norm(self.coords_f - self.coords_i)
        angle = np.arctan2(self.coords_f[1] - self.coords_i[1], self.coords_f[0] - self.coords_i[0])
        return length, angle

    def local_matrix (self):
        #L = self.length
        
        
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

        T = np.array([[c, s, 0, 0, 0, 0], 
                      [s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, s, c, 0],
                      [0, 0, 0, 0, 0, 1]])
        
        return np.dot(np.dot(T.T, k), T)
    
class StructureMatrix:

    def __init__ (self, elements,nodes, ndof=3):

        self.N_nodes = len(nodes)
        self.elements = elements
        self.N_elements = len(elements)
        self.ndof = ndof
        self.k_assembly = self.assembly()

    def o_matrix (self):
        #return np.zeros((self.N_nodes*self.ndof, self.N_nodes*self.ndof))
        #Aqui debo trabajar asi ya que estoy trabajando con algebre
        return sp.Matrix.zeros(self.N_nodes*self.ndof, self.N_nodes*self.ndof)

       
        
    def assembly (self):
        o_matrix = self.o_matrix()

        #Ahora debo sumar los elementos de la matriz de rigides, en base a las matrices globales de cada elemento
        for element in self.elements:
            k = element.k_global
         
            #Debo dividir la matriz en cuatro secciones iguales
            Q1 = k[:3, :3]  # Cuadrante superior izquierdo
            Q2 = k[:3, 3:]  # Cuadrante superior derecho
            Q3 = k[3:, :3]  # Cuadrante inferior izquierdo
            Q4 = k[3:, 3:]  # Cuadrante inferior derecho

            #Ahora debo sumar cada cuadrante a la martriz original
            
            n1 = element.n1
            n2 = element.n2

            # Sumar cada cuadrante a la matriz global en su ubicación correspondiente
            o_matrix[n1*3:(n1+1)*3, n1*3:(n1+1)*3] += Q1  # Cuadrante superior izquierdo en (n1, n1)
            o_matrix[n1*3:(n1+1)*3, n2*3:(n2+1)*3] += Q2  # Cuadrante superior derecho en (n1, n2)
            o_matrix[n2*3:(n2+1)*3, n1*3:(n1+1)*3] += Q3  # Cuadrante inferior izquierdo en (n2, n1)
            o_matrix[n2*3:(n2+1)*3, n2*3:(n2+1)*3] += Q4  # Cuadrante inferior derecho en (n2, n2)

        return o_matrix






        #Ahora debo sumar la matriz
        
        
        
    

    
#Defino algunos nodos
n1 = Node(0, np.array([0, 0]))
n2 = Node(1, np.array([3, 0]))
n3 = Node(2, np.array([0, 4]))

#Defino algunos elementos
e1 = Elements(n1, n2)


#Defino la matriz de estructura
sm = StructureMatrix([e1], [n1, n2])
sp.pprint(sm.k_assembly)