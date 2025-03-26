import numpy as np
from nodes import Node
from elements import Elements

class Structure:

    def __init__(self, name):
        self.name = name
        self.nodes = []
        self.elements = []
        self.espaciado_h = 9.14 #m
        self.Espaciado_v = [3.66, 
                    5.49, 
                    1.83, 
                    3.96 - 1.83, 
                    3.96, 
                    1.83, 
                    3.96 - 1.83, 
                    3.96, 
                    1.83, 
                    3.96 - 1.83, 
                    3.96, 
                    1.83, 
                    3.96 - 1.83, 
                    3.96]

#Defino los self.pisos en que hay cambios de seccion
        self.losas = [1, #Piso 1
                1, #Piso 2
                0, #Piso 3
                1, #Piso 4
                1, #Piso 5
                0, #Piso 6
                1, #Piso 7
                1, #Piso 8
                0, #Piso 9
                1, #Piso 10
                1, #Piso 11
                0, #Piso 12
                1, #Piso 13
                1] #Piso 14
        
#Ahora defino las secciones por piso [A, I]

#Seccion columnas (mm^2, mm^4)   I y D       Centrales   
        self.secciones_piso = [[[70322, 2264298957.44], [94839, 3417000000]], #Piso 1
                        [[70322, 2264298957.44], [94839, 3417000000]], #Piso 2
                        [[70322, 2264298957.44], [94839, 3417000000]], #Piso 3
                        [[70322, 2264298957.44], [86451, 2993000000]], #Piso 4
                        [[70322, 2264298957.44], [86451, 2993000000]], #Piso 5
                        [[70322, 2264298957.44], [86451, 2993000000]], #Piso 6
                        [[53742, 1598000000], [70322, 2264298957.44]], #Piso 7
                        [[53742, 1598000000], [70322, 2264298957.44]], #Piso 8
                        [[53742, 1598000000], [70322, 2264298957.44]], #Piso 9
                        [[48774, 1415000000], [53742, 1598000000]], #Piso 10
                        [[48774, 1415000000], [53742, 1598000000]], #Piso 11
                        [[48774, 1415000000], [53742, 1598000000]], #Piso 12
                        [[44193, 1252856592,26], [48774, 1415000000]], #Piso 13
                        [[44193, 1252856592,26], [48774, 1415000000]]  #Piso 14
                        ]

        
        #Defino los self.secciones_vigas en cada piso
        self.secciones_vigas = [[30323, 4062418717.759], #Piso 1
                [30323, 4062418717.759], #Piso 2
                0, #Piso 3
                [22064, 2052000000], #Piso 4
                [12968, 761703509.579], #Piso 5
                0, #Piso 6
                [22064, 2052000000], #Piso 7
                [12968, 761703509.579], #Piso 8
                0, #Piso 9
                [22064, 2052000000], #Piso 10
                [12968, 761703509.579], #Piso 11
                0, #Piso 12
                [22064, 2052000000], #Piso 13
                [12968, 761703509.579]] #Piso 14

        self.cargas_q = [26,
                    26,
                    0,
                    26,
                    26,
                    0,
                    26,
                    26,
                    0,
                    26,
                    26,
                    0,
                    26,
                    23.1]
        
        
        self.nodos_ancho = 6
        self.pisos = 14
        self.base_v = 0
        self.assembly()

    def assembly (self):
        for i in range(self.nodos_ancho):
            self.nodes.append(Node(i, np.array([i*self.espaciado_h, 0]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))


        #Defino los nodos de los self.pisos
        for j in range(1, self.pisos+1):
            espaciado_v = self.Espaciado_v[j-1]
            vertical = self.base_v + espaciado_v

            if j == 1:
                for i in range(self.nodos_ancho):
                    if i == 0:
                        self.nodes.append(Node(i+j*self.nodos_ancho, np.array([i*self.espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

                    elif i == self.nodos_ancho-1:
                        self.nodes.append(Node(i+j*self.nodos_ancho, np.array([i*self.espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

                    else:
                        self.nodes.append(Node(i+j*self.nodos_ancho, np.array([i*self.espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

            else:
                for i in range(self.nodos_ancho):

                    if i == 0:
                        self.nodes.append(Node(i+j*self.nodos_ancho, np.array([i*self.espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

                    else:
                        self.nodes.append(Node(i+j*self.nodos_ancho, np.array([i*self.espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

            #Ahora conecto verticalmente los self.pisos
            for i in range(self.nodos_ancho):
                if i == 0:
                    AI = self.secciones_piso[j-1][0]
                elif i == self.nodos_ancho-1:
                    AI = self.secciones_piso[j-1][0]
                else:
                    AI = self.secciones_piso[j-1][1]
                self.elements.append(Elements(self.nodes[i+(j-1)*self.nodos_ancho], self.nodes[i+j*self.nodos_ancho], AI=AI, dxdy=[0.8, 0]))

            if self.losas[j-1] == 1:
                AI = self.secciones_vigas[j-1]
                #Defino elementos horizontales que conectan los ultimos nodos creados
                for i in range(self.nodos_ancho-1):
                    nodos_actuales = len(self.nodes)
                    self.elements.append(Elements(self.nodes[nodos_actuales-(self.nodos_ancho) + i], self.nodes[nodos_actuales-self.nodos_ancho + i + 1 ], AI = AI, q=-1000*self.cargas_q[j-1], dxdy=[3, 0]))

            self.base_v = vertical