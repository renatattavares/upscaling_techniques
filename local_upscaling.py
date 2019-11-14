"""
Module of local upscaling technique in structured tridimensional meshes
"""

from import_impress import preprocessor
import numpy as np

class local_upscaling:
    def __init__(self):
        """
        Steps of local upscaling:

    -> Apply boundary conditions to local problems
    -> Assembly of local problems
    -> Solve local problems
    -> Calculate upscaled permeability or upscaled transmissibility
    -> Apply boundary conditions to global coarse problem
    -> Assembly of global coarse problem
    -> Solve global coarse problem
        """
        pass

    def set_boundary_conditions(self):
        pass

    def assembly(self):
        pass

    def solver(self):
        pass
