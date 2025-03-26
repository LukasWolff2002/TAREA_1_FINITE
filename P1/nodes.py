import numpy as np
class Node:

    def __init__ (self, n, coord, boundary, force_vector, vector_desplazamientos = [0,0,0], ndof=3):
        self.n = n
        self.coord = coord
        self.ndof = ndof
        self.idx = self.index()
        self.boundary = boundary
        self.force_vector = force_vector
        self.vector_desplazamientos = vector_desplazamientos
        self.def_vector = np.array([0.0,0.0,0.0])
        self._Wind(1)

    def index (self):
        #Debo generar un array desde n hasta n+ndof, con np.array
        return np.linspace((self.n*self.ndof), (self.n*self.ndof)+(self.ndof-1), 3, dtype=int)
 
    def _Wind (self, w):

        if self.coord[0] == 0:
            self.force_vector += np.array([w,0,0])
