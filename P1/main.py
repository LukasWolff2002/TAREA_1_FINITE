import numpy as np
from solver import Solver
from matrix_assembly import Assembly
from structure_assembly import Structure
from graph import plot_original_structure_all_forces, plot_deformed_structure, plotStructureWithMomentDiagram, plotStructureWithShearDiagram, plotStructureWithAxialDiagram
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

plotStructureWithMomentDiagram(structure_1.elements, escala_momentos=1e-7, show_structure=True)
plotStructureWithShearDiagram(structure_1.elements, escala_corte=1e-4, show_structure=False)
plotStructureWithAxialDiagram(structure_1.elements, escala_axial=1e-3, show_structure=True)