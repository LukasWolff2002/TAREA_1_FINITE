import numpy as np
from solver import Solver
from matrix_assembly import Assembly
from structure_assembly import Structure
from graph import plot_original_structure_all_forces, plot_deformed_structure
import matplotlib.pyplot as plt


#---------------------------------
#Estructura 1
#No tiene rigid offsets
#---------------------------------

structure_1 = Structure(dxdy = [0.0, 0.0])

M = Assembly(structure_1.nodes, structure_1.elements)

plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)

Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)
 
plot_deformed_structure(structure_1.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=10)

#---------------------------------
#Estructura 2
#Tiene rigid offsets
#---------------------------------

structure_2 = Structure(dxdy = [0.5, 0.0])

plot_original_structure_all_forces(structure_2.nodes, structure_2.elements)

M = Assembly(structure_2.nodes, structure_2.elements)

Solver(structure_2.nodes, structure_2.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

plot_deformed_structure(structure_2.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=10)
