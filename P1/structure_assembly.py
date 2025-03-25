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
        
#Ahora defino las secciones por piso

#Seccion (mm, mm)   I y D       Centrales   
        self.secciones_piso = [[[300, 300], [100, 200]], #Piso 1
                        [[300, 300], [100, 200]], #Piso 2
                        [[300, 300], [100, 200]], #Piso 3
                        [[300, 300], [100, 200]], #Piso 4
                        [[300, 300], [100, 200]], #Piso 5
                        [[300, 300], [100, 200]], #Piso 6
                        [[300, 300], [100, 200]], #Piso 7
                        [[300, 300], [100, 200]], #Piso 8
                        [[300, 300], [100, 200]], #Piso 9
                        [[300, 300], [100, 200]], #Piso 10
                        [[300, 300], [100, 200]], #Piso 11
                        [[300, 300], [100, 200]], #Piso 12
                        [[300, 300], [100, 200]], #Piso 13
                        [[300, 300], [100, 200]]  #Piso 14
                        ]

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
                    A = self.secciones_piso[j-1][0]
                elif i == self.nodos_ancho-1:
                    A = self.secciones_piso[j-1][0]
                else:
                    A = self.secciones_piso[j-1][1]
                self.elements.append(Elements(self.nodes[i+(j-1)*self.nodos_ancho], self.nodes[i+j*self.nodos_ancho], A=A, dxdy=[0.8, 0]))

            if self.losas[j-1] == 1:
                #Defino elementos horizontales que conectan los ultimos nodos creados
                for i in range(self.nodos_ancho-1):
                    nodos_actuales = len(self.nodes)
                    self.elements.append(Elements(self.nodes[nodos_actuales-(self.nodos_ancho) + i], self.nodes[nodos_actuales-self.nodos_ancho + i + 1 ], q=-1000*self.cargas_q[j-1], dxdy=[3, 0]))

            self.base_v = vertical