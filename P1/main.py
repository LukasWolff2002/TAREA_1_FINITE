import numpy as np
from solver import Solver
from matrix_assembly import Assembly
from structure_assembly import Structure
from graph import plot_original_structure_all_forces, plot_deformed_structure, plot_structure_with_local_displacements
import matplotlib.pyplot as plt

structure_1 = Structure("Estructura 1")

plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)

M = Assembly(structure_1.nodes, structure_1.elements)

Des = Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

#Ya,ahora debo obtener las distinstas cosas

#Desplazamientos

