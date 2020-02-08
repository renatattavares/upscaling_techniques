"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems

class LocalUpscaling:

    def __init__(self, mesh_file = None, boundary_condition_type = None):
        initial_time = time.time()

        print("\n########## Local upscaling class initialized ##########")
        lp = LocalProblems(mesh_file, boundary_condition_type)

        #self.number_volumes_local_problem = self.lp.number_volumes_local_problem
        # self.number_coarse_volumes = self.lp.number_coarse_volumes
        # self.x = self.lp.x
        # self.y = self.lp.y
        # self.z = self.lp.z
        #
        # # Parameters that need to be accessed in IMEX dataset
        # self.viscosity = 1
        # self.face_area = 1

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        for i in self.lp.direction_string:
            self.direction = self.lp.directions_dictionary.get(i)

            if np.array_equal(self.direction, self.x) is True:
                correct_volumes = self.lp.correct_volumes_x
            elif np.array_equal(self.direction, self.y) is True:
                correct_volumes = self.lp.correct_volumes_y
            elif np.array_equal(self.direction, self.z) is True:
                correct_volumes = self.lp.correct_volumes_z

            for j in range(self.number_coarse_volumes):

                # Check if faces are in the correct direction
                adjacent_faces = np.unique(self.coarse.elements[j].volumes.adjacencies[correct_volumes[j]].flatten())
                faces_normal = self.coarse.elements[j].faces.normal[adjacent_faces.flatten()]
                cross_product = np.cross(faces_normal, self.direction)
                norm_cross_product = np.linalg.norm(cross_product, axis = 1)
                parcial_correct_faces = np.isin(norm_cross_product, 0)
                index_correct_faces = np.where(parcial_correct_faces == True)[0]
                parcial_correct_faces = adjacent_faces[index_correct_faces]

                # Check if faces are internal
                internal_faces = self.coarse.elements[j].faces.internal
                correct_faces = np.isin(parcial_correct_faces, internal_faces)
                index_correct_faces = np.where(correct_faces == True)[0]
                correct_faces = parcial_correct_faces[index_correct_faces]

                equivalent_permeability = self.coarse.elements[j].equivalent_permeability[correct_faces]
                adjacent_volumes = self.coarse.elements[j].faces.bridge_adjacencies(correct_faces, 2, 3)
                global_id_adjacent_volumes = self.coarse.elements[j].volumes.father_id[adjacent_volumes.flatten()]
                pressures = self.mesh.pressure_x[global_id_adjacent_volumes]


                flow_rate = (-1/self.viscosity)*self.equivalent_permeability*self.face_area
                total_flow_rate = flow_rate.sum()
                #self.effective_permeability = total_flow_rate*/

    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.father_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.arange(len(porosity))
            volumes[:] = 1
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()
