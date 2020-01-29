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
        self.x = self.lp.x
        self.y = self.lp.y
        self.z = self.lp.z

        # Parameters that need to be accessed in IMEX dataset
        self.viscosity = 1

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        # for i in self.lp.direction_string:
        #     self.direction = self.lp.directions_dictionary.get(i)
        #
        #     if np.array_equal(self.direction, self.x) is True:
        #         correct_volumes = self.lp.correct_volumes_x
        #     elif np.array_equal(self.direction, self.y) is True:
        #         correct_volumes = self.lp.correct_volumes_y
        #     elif np.array_equal(self.direction, self.z) is True:
        #         correct_volumes = self.lp.correct_volumes_z

        correct_volumes = self.lp.correct_volumes_x
        self.direction = self.lp.x

        for j in range(1):
            # Check if faces are in the correct direction
            self.adjacent_faces = np.unique(self.coarse.elements[j].volumes.adjacencies[correct_volumes[j]].flatten())
            self.faces_normal = self.coarse.elements[j].faces.normal[self.adjacent_faces.flatten()]
            self.cross_product = np.cross(self.faces_normal, self.direction)
            self.norm_cross_product = np.linalg.norm(self.cross_product, axis = 1)
            # self.correct_faces = np.isin(self.norm_cross_product, 1)
            # self.index_correct_faces = np.where(self.correct_faces == True)[0]
            # self.correct_faces = self.adjacent_faces.flatten()[self.index_correct_faces]
            #
            # self.global_ids_correct_faces = self.coarse.elements[j].faces.father_id[self.correct_faces]


                # # Check if faces are internal
                # internal_faces = self.coarse.elements[j].faces.internal
                # correct_faces = np.isin(correct_faces, internal_faces)
                # index_correct_faces = np.where(correct_faces == False)[0]
                # correct_faces = adjacent_faces.flatten()[index_correct_faces]
                #
                # global_ids_correct_faces = self.coarse.elements[j].faces.father_id[correct_faces]
                # self.mesh.teste[global_ids_correct_faces] = 1000
                # #flow_rate = (-1/self.viscosity)*self.mesh.permeability[i]*
                self.mesh.teste[global_id] = 1000

    def upscale_porosity(self):
        pass

    def boundary_conditions(self):
        # Read data
        pass
