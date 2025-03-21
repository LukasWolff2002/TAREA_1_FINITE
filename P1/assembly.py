import numpy as np

class Assembly:

    def __init__ (self, nodes, elements, ndof=3):
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
            o_matrix[n1.n*3:(n1.n+1)*3, n1.n*3:(n1.n+1)*3] += Q1  # Cuadrante superior izquierdo en (n1, n1)
            o_matrix[n1.n*3:(n1.n+1)*3, n2.n*3:(n2.n+1)*3] += Q2  # Cuadrante superior derecho en (n1, n2)
            o_matrix[n2.n*3:(n2.n+1)*3, n1.n*3:(n1.n+1)*3] += Q3  # Cuadrante inferior izquierdo en (n2, n1)
            o_matrix[n2.n*3:(n2.n+1)*3, n2.n*3:(n2.n+1)*3] += Q4  # Cuadrante inferior derecho en (n2, n2)


        return o_matrix

    def extraer_submatriz(self, fila_cond, col_cond):
        o_matrix = self.k_assembly
        filas = set()
        cols = set()

        for node in self.nodes:
            for i, b in enumerate(node.boundary):
                if fila_cond(b):
                    filas.add(node.idx[i])
                if col_cond(b):
                    cols.add(node.idx[i])

        filas_restantes = [i for i in range(o_matrix.shape[0]) if i not in filas]
        cols_restantes = [j for j in range(o_matrix.shape[1]) if j not in cols]

        return o_matrix[np.ix_(filas_restantes, cols_restantes)]
    
    def final_matrix(self):
        return self.extraer_submatriz(lambda b: b == 1, lambda b: b == 1)

    def matriz_kfc(self):
        return self.extraer_submatriz(lambda b: b == 1, lambda b: b == 0)

    def matriz_kcf(self):
        return self.extraer_submatriz(lambda b: b == 0, lambda b: b == 1)

    def matriz_kcc(self):
        return self.extraer_submatriz(lambda b: b == 0, lambda b: b == 0)

