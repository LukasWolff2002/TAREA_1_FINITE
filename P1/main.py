import numpy as np
from nodes import Node
from elements import Elements
from desplazamientos import Desplazamientos
from assembly import Assembly

fc = 28
E = 4700 * np.sqrt(fc)
A = 300 * 300
I = 300 * 300 ** 3 / 12

A = A * 1000

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
sm = Assembly([e1, e2, e3], [n1, n2, n3, n4])

Des = Desplazamientos([n1, n2, n3, n4], sm.kff_matrix, sm.kfc_matrix, sm.kcf_matrix, sm.kcc_matrix)

print(Des.uf_v)

import matplotlib.pyplot as plt

def plot_structure(nodes, elements, desplazamientos, scale=1000):



