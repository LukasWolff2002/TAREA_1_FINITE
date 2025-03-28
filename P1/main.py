import numpy as np
from solver import Solver
from matrix_assembly import Assembly
from structure_assembly import Structure
from graph import plot_original_structure_all_forces, plot_deformed_structures, plotStructureWithMomentDiagram, plotStructureWithShearDiagram, plotStructureWithAxialDiagram
import matplotlib.pyplot as plt


#---------------------------------
#Sin offsets
#---------------------------------

structure_1 = Structure(nombre='Caso A, sin offsets rigidos',dxdy = [0, 0.0], caso='a')
structure_2 = Structure(nombre='Caso B, sin offsets rigidos',dxdy = [0, 0.0], caso='b')
structure_3 = Structure(nombre='Caso C, sin offsets rigidos',dxdy = [0, 0.0], caso='c')


M = Assembly(structure_1.nodes, structure_1.elements)
Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_2.nodes, structure_2.elements)
Solver(structure_2.nodes, structure_2.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_3.nodes, structure_3.elements)
Solver(structure_3.nodes, structure_3.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

plot_deformed_structures([structure_1, structure_2, structure_3], text=False, nodes=True, nodes_labels=False, deformada=True, Escala=[750, 100, 100])

plotStructureWithMomentDiagram([structure_1, structure_2, structure_3], Escala=[1e-5, 1e-6, 1e-6], show_structure=True)
plotStructureWithShearDiagram([structure_1, structure_2, structure_3], Escala=[1e-5, 1e-6, 1e-6], show_structure=True)
plotStructureWithAxialDiagram([structure_1, structure_2, structure_3], Escala=[1e-6, 1e-6, 1e-6], show_structure=True)

#---------------------------------
#Con offsets
#---------------------------------
structure_1 = Structure(nombre='Caso A, con offsets rigidos', dxdy = [0.5, 0.0], caso='a')
structure_2 = Structure(nombre='Caso B, con offsets rigidos', dxdy = [0.5, 0.0], caso='b')
structure_3 = Structure(nombre='Caso C, con offsets rigidos', dxdy = [0.5, 0.0], caso='c')

M = Assembly(structure_1.nodes, structure_1.elements)
Solver(structure_1.nodes, structure_1.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_2.nodes, structure_2.elements)
Solver(structure_2.nodes, structure_2.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

M = Assembly(structure_3.nodes, structure_3.elements)
Solver(structure_3.nodes, structure_3.elements, M.kff_matrix, M.kfc_matrix, M.kcf_matrix, M.kcc_matrix)

plot_deformed_structures([structure_1, structure_2, structure_3], text=False, nodes=True, nodes_labels=False, deformada=True, Escala=[750, 500, 500])

plotStructureWithMomentDiagram([structure_1, structure_2, structure_3], Escala=[1e-5, 1e-6, 1e-6], show_structure=True)
plotStructureWithShearDiagram([structure_1, structure_2, structure_3], Escala=[1e-6, 1e-6, 1e-6], show_structure=True)
plotStructureWithAxialDiagram([structure_1, structure_2, structure_3], Escala=[1e-6, 1e-6, 1e-6], show_structure=True)

