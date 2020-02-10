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
        self.px = lp.pressure_x
        self.py = lp.pressure_y
        self.pz = lp.pressure_z
        self.mesh = lp.mesh
        self.coarse = lp.coarse
        self.number_coarse_volumes = lp.number_coarse_volumes
        self.cx = lp.correct_volumes_x
        self.cy = lp.correct_volumes_y
        self.cz = lp.correct_volumes_z
        self.lp = lp

        self.effective_permeability = np.array([])
        self.upscale_permeability()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        for i in self.lp.direction_string:

            self.direction = self.lp.directions_dictionary.get(i)

            if np.array_equal(self.direction, self.lp.x) is True:
                correct_volumes = self.cx
                direction_number = 0
            elif np.array_equal(self.direction, self.lp.y) is True:
                correct_volumes = self.cy
                direction_number = 1
            elif np.array_equal(self.direction, self.lp.z) is True:
                correct_volumes = self.cz
                direction_number = 2

            for j in range(self.number_coarse_volumes):
                # Check if faces are in correct direction
                global_id_faces = self.coarse.elements[j].faces.global_id[:]
                local_id_faces = self.coarse.elements[j].faces.all
                parallel_direction = self.mesh.parallel_direction[global_id_faces]
                correct_faces = np.isin(parallel_direction, direction_number)
                index_correct_faces = np.where(correct_faces == True)[0]
                partial_correct_faces = global_id_faces[index_correct_faces]
                partial_correct_faces_local = local_id_faces[index_correct_faces]

                # Check if faces are internal
                internal_faces = self.coarse.elements[j].faces.internal
                global_id_faces = self.coarse.elements[j].faces.global_id[internal_faces]
                correct_faces = np.isin(partial_correct_faces, global_id_faces)
                index_correct_faces = np.where(correct_faces == True)[0]
                correct_faces = partial_correct_faces[index_correct_faces]
                correct_faces_local = partial_correct_faces_local[index_correct_faces]

                equivalent_permeability = self.mesh.equivalent_permeability[correct_faces]
                adjacent_volumes = self.coarse.elements[j].faces.bridge_adjacencies(correct_faces_local, 2, 3)
                global_id_adjacent_volumes = self.coarse.elements[j].volumes.father_id[adjacent_volumes.flatten()]
                pressures = self.mesh.pressure_x[global_id_adjacent_volumes]


                flow_rate = (-1/self.viscosity)*self.equivalent_permeability*1*1
                total_flow_rate = flow_rate.sum()
                print(total_flow_rate)
                self.effective_permeability = np.append(self.effective_permeability, total_flow_rate/(1*1))

    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.global_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.ones(len(porosity))
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()
