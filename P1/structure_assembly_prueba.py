import numpy as np
from nodes import Node
from elements import Elements

class Structure:

    def __init__(self, name):
        self.name = name
        self.nodes = []
        self.elements = []
        self.espaciado_h = 9.14  # m
        self.Espaciado_v = [
            3.66, 5.49, 1.83, 3.96 - 1.83, 3.96, 1.83, 
            3.96 - 1.83, 3.96, 1.83, 3.96 - 1.83, 3.96, 
            1.83, 3.96 - 1.83, 3.96
        ]

        # Losas por piso
        self.losas = [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]

        # Secciones de las columnas (6 tipos diferentes) [A, I]
        self.secciones_columnas = {
            1: [70322, 2264298957.44],  # Sección 1
            2: [94839, 3417000000],     # Sección 2
            3: [86451, 2993000000],     # Sección 3
            4: [53742, 1598000000],     # Sección 4
            5: [48774, 1415000000],     # Sección 5
            6: [44193, 1252856592.26]   # Sección 6
        }
        
        # Secciones de las vigas (según la imagen proporcionada)
        self.secciones_vigas = [
            [30323, 4062418717.759],  # W36x160
            [22387, 2456000000],      # W33x118
            [22064, 2052000000],      # W30x116
            [12968, 761703509.579],   # W24x68
        ]

        self.cargas_q = [26, 26, 0, 26, 26, 0, 26, 26, 0, 26, 26, 0, 26, 23.1]
        
        self.nodos_ancho = 6
        self.pisos = 14
        self.base_v = 0
        self.assembly()

    def assembly(self):
        """
        Ensambla la estructura, creando nodos y elementos (columnas y vigas) con las propiedades
        de sección correspondientes.
        """
        # Crear los nodos base
        for i in range(self.nodos_ancho):
            self.nodes.append(Node(i, np.array([i * self.espaciado_h, 0]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

        # Crear los nodos de los pisos
        for j in range(1, self.pisos + 1):
            espaciado_v = self.Espaciado_v[j - 1]
            vertical = self.base_v + espaciado_v

            # Crear nodos en el piso
            for i in range(self.nodos_ancho):
                if i == 0 or i == self.nodos_ancho - 1:
                    self.nodes.append(Node(i + j * self.nodos_ancho, np.array([i * self.espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))
                else:
                    self.nodes.append(Node(i + j * self.nodos_ancho, np.array([i * self.espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

            # Conectar verticalmente los pisos
            for i in range(self.nodos_ancho):
                tipo_seccion = i % 6 + 1  # Esto asigna 1, 2, 3, ..., 6
                A = self.secciones_columnas[tipo_seccion]  # Sección de la columna para el tipo de nodo
                self.elements.append(Elements(self.nodes[i + (j - 1) * self.nodos_ancho], self.nodes[i + j * self.nodos_ancho], AI=A, dxdy=[0.8, 0]))

            # Agregar vigas horizontales si es necesario
            if self.losas[j - 1] == 1:
                for i in range(self.nodos_ancho - 1):
                    nodos_actuales = len(self.nodes)
                    # Asignar la sección de viga correspondiente
                    viga_seccion = self.secciones_vigas[i % len(self.secciones_vigas)]
                    self.elements.append(Elements(self.nodes[nodos_actuales - (self.nodos_ancho) + i], self.nodes[nodos_actuales - self.nodos_ancho + i + 1], q=-1000 * self.cargas_q[j - 1], dxdy=[3, 0]))

            # Actualizar la base vertical para el siguiente piso
            self.base_v = vertical
