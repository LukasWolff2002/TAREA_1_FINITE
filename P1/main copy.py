import numpy as np
from nodes import Node
from elements import Elements
from solver import Solver
from matrix_assembly import Assembly
from graph import plot_original_structure_all_forces, plot_deformed_structure, plot_structure_with_local_displacements
import matplotlib.pyplot as plt

nodes = []
elements = []

espaciado_h = 9.14 #m
Espaciado_v = [3.66, 
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

#Defino los pisos en que hay cambios de seccion
losas = [1, #Piso 1
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
secciones_piso = [[[300, 300], [100, 200]], #Piso 1
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

cargas_q = [26,
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
            
            

nodos_ancho = 6
pisos = 14

base_v = 0
#Defino los nodos base
for i in range(nodos_ancho):
    nodes.append(Node(i, np.array([i*espaciado_h, 0]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))


#Defino los nodos de los pisos
for j in range(1, pisos+1):
    espaciado_v = Espaciado_v[j-1]
    vertical = base_v + espaciado_v

    if j == 1:
        for i in range(nodos_ancho):
            if i == 0:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

            elif i == nodos_ancho-1:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([1, 1, 0]), np.array([0.0, 0.0, 0.0])))

            else:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

    else:
        for i in range(nodos_ancho):

            if i == 0:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

            else:
                nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

    #Ahora conecto verticalmente los pisos
    for i in range(nodos_ancho):
        if i == 0:
            A = secciones_piso[j-1][0]
        elif i == nodos_ancho-1:
            A = secciones_piso[j-1][0]
        else:
            A = secciones_piso[j-1][1]
        elements.append(Elements(nodes[i+(j-1)*nodos_ancho], nodes[i+j*nodos_ancho], A=A))

    if losas[j-1] == 1:
        #Defino elementos horizontales que conectan los ultimos nodos creados
        for i in range(nodos_ancho-1):
            nodos_actuales = len(nodes)
            elements.append(Elements(nodes[nodos_actuales-(nodos_ancho) + i], nodes[nodos_actuales-nodos_ancho + i + 1 ], q=-1000*cargas_q[j-1]))

    base_v = vertical

#lot_original_structure_all_forces(nodes, elements)

M = Assembly(nodes, elements)

Des = Solver(nodes, elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

#Ya,ahora debo obtener las distinstas cosas

#Desplazamientos

for node in nodes:
    print(f"Nodo {node.n}")
    print(f"Desplazamientos: {node.def_vector}")

def extractDisplacements(self): #Extraer los desplazamientos de los nodos
        Tgo = self.tgo
        un1 = np.array(self.n1.def_vector)
        un2 = np.array(self.n2.def_vector)
        u_global = np.concatenate([un1, un2])
        Tlg = self.localGlobalTransformation()
        u_local = Tlg @ u_global
        Tbl = self.basicLocalTransformation()
        u_basic = Tbl @ u_local
        return u_global, u_local, u_basic

for element in elements:
    element.ug, element.ul, element.ub = extractDisplacements(element)
    #print(f"Elemento {element.n1.n} - {element.n2.n}")
    #print(f"Desplazamientos globales: {element.ug}") 
    #print(f"Desplazamientos locales: {element.ul}")
    #print(f"Desplazamientos basicos: {element.ub}")

def calculateBasicForces(self):
        u_basic = self.ub
        Kb = self.kb
        f_basic = Kb @ u_basic
        return f_basic


def calculateLocalForces(self):
    u_local = self.ul
    Kl = self.kl
    f_local = Kl @ u_local
    return f_local


def forceRecovery(self):
        f_basic = calculateBasicForces(self)
        f_local = calculateLocalForces(self)
        Tlg = self.localGlobalTransformation()
        f_global = Tlg.T @ f_local

        return f_basic, f_local, f_global

for element in elements:
    element.fb, element.fl, element.fg = forceRecovery(element)

    print(f"Elemento {element.n1.n} - {element.n2.n}")
    print(f"Fuerzas basicas: {element.fb}")
    print(f"Fuerzas locales: {element.fl}")
    print(f"Fuerzas globales: {element.fg}")