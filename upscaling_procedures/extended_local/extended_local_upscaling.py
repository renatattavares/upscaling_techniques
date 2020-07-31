"""
Module of extended local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.extended_local.extended_local_problems import ExtendedLocalProblems

class ExtendedLocalUpscaling(ExtendedLocalProblems):

    def __init__(self, mesh_file = None, dataset = None):

        initial_time = time.time()

        print('\n##### ExtendedLocalUpscaling class initialized #####')

        super().__init__(mesh_file, dataset)
        self.center_distance_walls()
        self.areas()
        #self.upscale_permeability_porosity()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def center_distance_walls(self):

        self.center_distance_walls = []

        for cv in range(self.number_coarse_volumes):
            volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(cv)
            distance = []
            min_coord = self.mesh.volumes.center[volumes_global_ids].min(axis = 0)
            max_coord = self.mesh.volumes.center[volumes_global_ids].max(axis = 0)
            distance.append(max_coord[0] - min_coord[0])
            distance.append(max_coord[1] - min_coord[1])
            distance.append(max_coord[2] - min_coord[2])
            self.center_distance_walls.append(distance)
