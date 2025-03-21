import numpy as np
from nodes import Node
from elements import Elements
from desplazamientos import Desplazamientos
from assembly import Assembly
from graph import plot_original_structure_all_forces, plot_deformed_structure
import matplotlib.pyplot as plt

nodes = []
elements = []

espaciado_h = 7 #m
Espaciado_v = [3.66, 
               5.49, 
               (3.96)/2, 
               (3.96)/2, 
               3.96, 
               (3.96)/2, 
               (3.96/2), 
               3.96, 
               (3.96)/2, 
               (3.96/2), 
               3.96, 
               (3.96)/2, 
               (3.96/2),
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

#Seccion (mm, mm)   Izquierda   Centrales    Derecha
secciones_piso = [[[300, 300], [100, 100], [100, 100]], #Piso 1
                  [[300, 300], [100, 100], [100, 100]], #Piso 2
                  [[300, 300], [100, 100], [100, 100]], #Piso 3
                  [[300, 300], [100, 100], [100, 100]], #Piso 4
                  [[300, 300], [100, 100], [100, 100]], #Piso 5
                  [[300, 300], [100, 100], [100, 100]], #Piso 6
                  [[300, 300], [100, 100], [100, 100]], #Piso 7
                  [[300, 300], [100, 100], [100, 100]], #Piso 8
                  [[300, 300], [100, 100], [100, 100]], #Piso 9
                  [[300, 300], [100, 100], [100, 100]], #Piso 10
                  [[300, 300], [100, 100], [100, 100]], #Piso 11
                  [[300, 300], [100, 100], [100, 100]], #Piso 12
                  [[300, 300], [100, 100], [100, 100]], #Piso 13
                  [[300, 300], [100, 100], [100, 100]]  #Piso 14
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
            
            

nodos_ancho = 3 #6
pisos = 14

base_v = 0
#Defino los nodos base
for i in range(nodos_ancho):
    nodes.append(Node(i, np.array([i*espaciado_h, 0]), np.array([1, 1, 1]), np.array([0.0, 0.0, 0.0])))


#Defino los nodos de los pisos
for j in range(1, pisos+1):
    espaciado_v = Espaciado_v[j-1]
    vertical = base_v + espaciado_v

   
    for i in range(nodos_ancho):
        nodes.append(Node(i+j*nodos_ancho, np.array([i*espaciado_h, vertical]), np.array([0, 0, 0]), np.array([0.0, 0.0, 0.0])))

    #Ahora conecto verticalmente los pisos
    for i in range(nodos_ancho):
        if i == 0:
            A = secciones_piso[j-1][0]
        elif i == nodos_ancho-1:
            A = secciones_piso[j-1][2]
        else:
            A = secciones_piso[j-1][1]
        elements.append(Elements(nodes[i+(j-1)*nodos_ancho], nodes[i+j*nodos_ancho], A=A))

    if losas[j-1] == 1:
        #Defino elementos horizontales que conectan los ultimos nodos creados
        for i in range(nodos_ancho-1):
            nodos_actuales = len(nodes)
            elements.append(Elements(nodes[nodos_actuales-(nodos_ancho) + i], nodes[nodos_actuales-nodos_ancho + i + 1 ], q=-1000*cargas_q[j-1]))

    base_v = vertical

#Ahora tomo la masa total de la estructra

M = Assembly(nodes, elements)

Des = Desplazamientos(nodes, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

# Ejecutar la función para graficar la estructura c on todas las fuerzas aplicadas
plot_original_structure_all_forces(nodes, elements)



# Ejecutar la función para graficar la estructura desplazada con escala ajustable
plot_deformed_structure(nodes, elements, scale=100)

for node in nodes:
    print(f"Nodo {node.n} - Desplazamientos: {node.def_vector}")
