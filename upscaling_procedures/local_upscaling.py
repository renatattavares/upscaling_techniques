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
        print("\n########## Local upscaling class initialized ##########")
        self.lp = local_problems_object(preprocessor, coarse_config, mesh_file, boundary_condition_type)
        self.mesh = self.lp.mesh
        self.coarse = self.lp.coarse
        self.number_volumes_local_problem = self.lp.number_volumes_local_problem
        self.number_coarse_volumes = self.lp.number_coarse_volumes
        #self.permeability = self.lp.permeability

    def upscaled_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """



    def assembly(self):
        pass

    def boundary_conditions(self):
        # Read data
        pass

    def solver(self):
        pass
