import numpy as np

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
        self.rc_v = self.rc_vector()

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
