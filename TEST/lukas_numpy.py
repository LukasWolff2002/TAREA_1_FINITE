import numpy as np
import sympy as sp


fc = 28
E = 4700 * np.sqrt(fc)
A = 300 * 300
I = 300 * 300 ** 3 / 12

A = A * 1

class Node:

    def __init__ (self, n, coord, boundary, force_vector, vector_desplazamientos = [0,0,0], ndof=3):
        self.n = n
        self.coord = coord
        self.ndof = ndof
        self.idx = self.index()
        self.boundary = boundary
        self.force_vector = force_vector
        self.vector_desplazamientos = vector_desplazamientos

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
    
class StructureMatrix:

    def __init__ (self, elements,nodes, ndof=3):
        self.nodes = nodes
        self.N_nodes = len(nodes)
        self.elements = elements
        self.N_elements = len(elements)
        self.ndof = ndof
        self.k_assembly = self.assembly()
        self.kff_matrix = self.final_matrix()
        self.kfc_matrix = self.matriz_kfc()
        self.kcf_matrix = self.matriz_kcf()
        self.kcc_matrix = self.matriz_kcc()

    def o_matrix (self):
        return np.zeros((self.N_nodes*self.ndof, self.N_nodes*self.ndof))
        
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

    def final_matrix(self):
        o_matrix = self.k_assembly  # Suponemos que es una matriz NumPy

        # Diccionarios para almacenar los índices a eliminar
        filas_a_eliminar = set()
        columnas_a_eliminar = set()

        # Determinar qué filas y columnas deben eliminarse
        for nodes in self.nodes:
            for i, b in enumerate(nodes.boundary):
                if b == 1:
                    filas_a_eliminar.add(nodes.idx[i])
                    columnas_a_eliminar.add(nodes.idx[i])

        # Crear listas con los índices que queremos mantener
        filas_restantes = [i for i in range(o_matrix.shape[0]) if i not in filas_a_eliminar]
        columnas_restantes = [j for j in range(o_matrix.shape[1]) if j not in columnas_a_eliminar]

        # Extraer la matriz sin esas filas y columnas
        o_matrix = o_matrix[np.ix_(filas_restantes, columnas_restantes)]

        return o_matrix

    
  

    def matriz_kfc(self):
        o_matrix = self.k_assembly  # Suponemos que es una matriz NumPy

        filas_a_eliminar = set()
        columnas_a_eliminar = set()

        # Determinar qué filas y columnas deben eliminarse
        for nodes in self.nodes:
            for i, b in enumerate(nodes.boundary):
                if b == 1:
                    filas_a_eliminar.add(nodes.idx[i])  # Eliminar columna si boundary es 1
                if b == 0:
                    columnas_a_eliminar.add(nodes.idx[i])  # Eliminar fila si boundary es 0


        # Crear listas con los índices que queremos mantener
        filas_restantes = [i for i in range(o_matrix.shape[0]) if i not in filas_a_eliminar]
        columnas_restantes = [j for j in range(o_matrix.shape[1]) if j not in columnas_a_eliminar]

        # Extraer la submatriz con las filas y columnas restantes
        o_matrix = o_matrix[np.ix_(filas_restantes, columnas_restantes)]

        return o_matrix
    
    def matriz_kcf(self):
        o_matrix = self.k_assembly  # Suponemos que es una matriz NumPy

        filas_a_eliminar = set()
        columnas_a_eliminar = set()

        # Determinar qué filas y columnas deben eliminarse
        for nodes in self.nodes:
            for i, b in enumerate(nodes.boundary):
                if b == 1:
                    columnas_a_eliminar.add(nodes.idx[i])  # Eliminar columna si boundary es 1
                if b == 0:
                    filas_a_eliminar.add(nodes.idx[i])  # Eliminar fila si boundary es 0

        # Crear listas con los índices que queremos mantener
        filas_restantes = [i for i in range(o_matrix.shape[0]) if i not in filas_a_eliminar]
        columnas_restantes = [j for j in range(o_matrix.shape[1]) if j not in columnas_a_eliminar]

        # Extraer la submatriz con las filas y columnas restantes
        o_matrix = o_matrix[np.ix_(filas_restantes, columnas_restantes)]

        return o_matrix
    
    def matriz_kcc(self):
        o_matrix = self.k_assembly  # Suponemos que es una matriz NumPy

        # Diccionarios para almacenar los índices a eliminar
        filas_a_eliminar = set()
        columnas_a_eliminar = set()

        # Determinar qué filas y columnas deben eliminarse
        for nodes in self.nodes:
            for i, b in enumerate(nodes.boundary):
                if b == 0:
                    filas_a_eliminar.add(nodes.idx[i])
                    columnas_a_eliminar.add(nodes.idx[i])

        # Crear listas con los índices que queremos mantener
        filas_restantes = [i for i in range(o_matrix.shape[0]) if i not in filas_a_eliminar]
        columnas_restantes = [j for j in range(o_matrix.shape[1]) if j not in columnas_a_eliminar]

        # Extraer la matriz sin esas filas y columnas
        o_matrix = o_matrix[np.ix_(filas_restantes, columnas_restantes)]

        return o_matrix


class Desplazamientos:

    def __init__(self, nodes, kff, kfc, kcf, kcc):
        self.nodes = nodes
        self.kff = kff
        self.kfc = kfc
        self.kcf = kcf
        self.kcc = kcc
        self.ff = self.ff_vector()
        self.uc = self.uc_vector()
        self.uf_v = self.uf_vector()
        self.rc = self.rc_vector()

    def ff_vector(self):
        #Debo generar un vector de fuerzas

        ff = []
        
        for nodes in self.nodes:
            for i, b in enumerate(nodes.boundary):
                if b == 1:
                    pass
                else:
                    ff.append(nodes.force_vector[i])

        ff = np.array(ff).reshape(-1, 1)
        return ff
    
    def uc_vector(self):
        #Debo generar un vector de desplazamientos

        uc = []
        
        for nodes in self.nodes:
            for i, b in enumerate(nodes.boundary):
                if b == 1:
                    uc.append(0)
                else:
                    pass

        uc = np.array(uc).reshape(-1, 1)

        return uc
    
    
    def uf_vector (self):
        uf = np.linalg.inv(self.kff) @ (self.ff - self.kfc @ self.uc)
        return uf
    
    def rc_vector (self):
        rc = self.kcf @ self.uf_v + self.kcc @ self.uc
        return rc

#Defino algunos nodos
n1 = Node(0, np.array([0, 0]), np.array([1, 1 ,1]), np.array([0, 0, 0]))
n2 = Node(1, np.array([2, 3]), np.array([0, 0, 0]), np.array([29419.95, -9806.65, 0]))
n3 = Node(2, np.array([6, 3]), np.array([0, 0, 0]), np.array([0, -9806.65, 0]))
n4 = Node(3, np.array([6, 0]), np.array([1, 1, 0]), np.array([0, 0, 0]))


#Defino algunos elementos
e1 = Elements(n1, n2)
e2 = Elements(n2, n3)
e3 = Elements(n3, n4) 


#Defino la matriz de estructura
sm = StructureMatrix([e1, e2, e3], [n1, n2, n3, n4])

Des = Desplazamientos([n1, n2, n3, n4], sm.kff_matrix, sm.kfc_matrix, sm.kcf_matrix, sm.kcc_matrix)

print(Des.uf_v)

print(Des.rc)








