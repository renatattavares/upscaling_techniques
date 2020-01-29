"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from local_problems.local_problems import LocalProblems

class LocalUpscaling:

    def __init__(self, mesh_file = None, boundary_condition_type = None):
        initial_time = time.time()

        print("\n########## Local upscaling class initialized ##########")
        self.lp = LocalProblems(mesh_file, boundary_condition_type)
        self.mesh = self.lp.mesh
        self.coarse = self.lp.coarse
        self.number_volumes_local_problem = self.lp.number_volumes_local_problem
        self.number_coarse_volumes = self.lp.number_coarse_volumes

        # Parameters that need to be accessed
        self.viscosity = 1



        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        #flow_rate = (-1/self.viscosity)*self.mesh.permeability[i]*
        pass

    def upscale_porosity(self):
        pass

    def boundary_conditions(self):
        # Read data
        pass
