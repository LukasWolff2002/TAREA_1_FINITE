import numpy as np
from solver import Solver
from matrix_assembly import Assembly
from structure_assembly import Structure
from graph import plot_original_structure_all_forces, plot_deformed_structure, plot_deformed_structures, plotStructureWithMomentDiagram, plotStructureWithShearDiagram, plotStructureWithAxialDiagram
import matplotlib.pyplot as plt


#---------------------------------
#Sin offsets
#---------------------------------

structure_1 = Structure(dxdy = [0, 0.0], caso='a')
structure_2 = Structure(dxdy = [0, 0.0], caso='b')
structure_3 = Structure(dxdy = [0, 0.0], caso='c')


M = Assembly(structure_1.nodes, structure_1.elements)
Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_2.nodes, structure_2.elements)
Solver(structure_2.nodes, structure_2.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_3.nodes, structure_3.elements)
Solver(structure_3.nodes, structure_3.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)
plot_original_structure_all_forces(structure_2.nodes, structure_2.elements)
plot_original_structure_all_forces(structure_3.nodes, structure_3.elements)

#plot_deformed_structure(structure_1.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=4000)
#plot_deformed_structure(structure_2.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=2000)
#plot_deformed_structure(structure_3.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=2000)

#plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)

plot_deformed_structures([structure_1, structure_2, structure_3], text=False, nodes=True, nodes_labels=False, deformada=True, escala=500)

plotStructureWithMomentDiagram(structure_1.elements, escala_momentos=1e-5, show_structure=True)
plotStructureWithShearDiagram(structure_1.elements, escala_corte=1e-5, show_structure=True)
plotStructureWithAxialDiagram(structure_1.elements, escala_axial=1e-6, show_structure=True)

#---------------------------------
#Con offsets
#---------------------------------
structure_1 = Structure(dxdy = [0.5, 0.0], caso='a')
structure_2 = Structure(dxdy = [0.5, 0.0], caso='b')
structure_3 = Structure(dxdy = [0.5, 0.0], caso='c')


M = Assembly(structure_1.nodes, structure_1.elements)
Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_2.nodes, structure_2.elements)
Solver(structure_2.nodes, structure_2.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_3.nodes, structure_3.elements)
Solver(structure_3.nodes, structure_3.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)
plot_original_structure_all_forces(structure_2.nodes, structure_2.elements)
plot_original_structure_all_forces(structure_3.nodes, structure_3.elements)

#plot_deformed_structure(structure_1.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=4000)
#plot_deformed_structure(structure_2.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=2000)
#plot_deformed_structure(structure_3.elements, text=False, nodes=True, nodes_labels=False, deformada=True, escala=2000)

#plot_original_structure_all_forces(structure_1.nodes, structure_1.elements)

plot_deformed_structures([structure_1, structure_2, structure_3], text=False, nodes=True, nodes_labels=False, deformada=True, escala=500)

plotStructureWithMomentDiagram(structure_1.elements, escala_momentos=1e-5, show_structure=True)
plotStructureWithShearDiagram(structure_1.elements, escala_corte=1e-5, show_structure=True)
plotStructureWithAxialDiagram(structure_1.elements, escala_axial=1e-6, show_structure=True)

