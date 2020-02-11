"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems

class LocalUpscaling(LocalProblems):

    def __init__(self, mesh_file = None, boundary_condition_type = None):
        initial_time = time.time()

        lp = LocalProblems(mesh_file, boundary_condition_type)

        #self.upscale_permeability()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """

    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.global_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.ones(len(porosity))
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()
