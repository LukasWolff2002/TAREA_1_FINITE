import numpy as np

class Solver:

    def __init__(self, nodes, elements, kff, kfc, kcf, kcc):
        self.nodes = nodes
        self.elements = elements
        self.kff = kff
        self.kfc = kfc
        self.kcf = kcf
        self.kcc = kcc
        self.ff = self.ff_vector()
        self.uc = self.uc_vector()
        self.uf_v = self.uf_vector()
        self.rc_v = self.rc_vector()
        self.node_elements_def_vector()


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

        uf_v = np.linalg.inv(self.kff) @ (self.ff - self.kfc @ self.uc)
        return uf_v
    
    def rc_vector (self):

        rc_v = self.kcf @ self.uf_v + self.kcc @ self.uc
        return rc_v
    
    def node_elements_def_vector (self):
        contador = 0  # índice para self.uf_v

        uf_v = self.uf_v

        for node in self.nodes:
            for i, b in enumerate(node.boundary):
                if b == 0:  # Solo si el GDL está libre
                    if i == 0 or i == 1:
                        node.def_vector[i] += uf_v[contador] /1000
                    else:
                        node.def_vector[i] += uf_v[contador]  
                    contador += 1

        for element in self.elements:
            #Con eso se cuales nodos conectan el elemento
            #Por lo tanto, puedo ensamblar el vector del elemento

            u = np.zeros(6)

            u[0:3] = element.n1.def_vector
            u[3:6] = element.n2.def_vector

            element.extractDisplacements(u)
            element.extractForces()

    


      


                    

