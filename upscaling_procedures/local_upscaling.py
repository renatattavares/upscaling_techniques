"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np

class LocalUpscaling:

    def __init__(self, local_problems_object, preprocessor, coarse_config, mesh_file = None, boundary_condition_type = None):
        """
            Steps of local upscaling:

        -> Assembly of local problems OK
        -> Apply boundary conditions to local problems OK
        -> Solve local problems OK
        -> Calculate upscaled permeability or upscaled transmissibility
        -> Assembly of global coarse problem
        -> Apply boundary conditions to global coarse problem
        -> Solve global coarse problem
        """
        lp = local_problems_object(preprocessor, coarse_config, mesh_file, boundary_condition_type)

    def upscaled_permeability(self):
        pass

    def upscaled_transmissibility(self):
        pass

    def assembly(self):
        pass

    def boundary_conditions(self):
        # Read data
        pass

    def solver(self):
        pass
