"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems

class LocalUpscaling(LocalProblems):

    def __init__(self, mesh_file = None, dataset = None):
        initial_time = time.time()

        super().__init__(mesh_file, dataset)
        self.center_distance_walls()
        self.areas()
        self.upscale_permeability_porosity()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability_porosity(self):
        """
        It calculates the effective permeability of a coarse volume
        """

        self.effective_permeability = []
        self.effective_porosity = []
        volume = self.length_elements[0]*self.length_elements[1]*self.length_elements[2]

        for cv in range(self.number_coarse_volumes):
            print('\nUpscaling of local problem {}'.format(cv))
            self.coarse_volume = cv
            general_transmissibility = self.assembly_local_problem()
            effective_permeability = []
            total_volume = volume*self.number_volumes_local_problem
            global_ids_volumes = self.coarse.elements[cv].volumes.father_id[:]
            porosity = self.mesh.porosity[global_ids_volumes]
            effective_porosity = (volume*porosity).sum()/total_volume

            for direction in self.direction_string:
                #print('in {} direction'.format(direction))
                self.direction = direction
                transmissibility, source, local_wall = self.set_boundary_conditions(general_transmissibility)
                pressures = self.solver(transmissibility, source)

                direction_number = self.directions_numbers.get(direction)
                center_distance_walls = self.center_distance_walls[cv][direction_number]
                area = self.areas[direction_number]
                global_wall = self.coarse.elements[cv].volumes.global_id[local_wall]
                local_adj = self.identify_adjacent_volumes_to_wall(cv, local_wall)
                global_adj = self.coarse.elements[cv].volumes.global_id[local_adj]
                center_wall = self.coarse.elements[cv].volumes.center[local_wall]
                center_adj = self.coarse.elements[cv].volumes.center[local_adj]
                pressure_wall = pressures[local_wall]
                pressure_adj = pressures[local_adj]
                pw = self.mesh.permeability[global_wall]
                permeability_wall = pw[:, direction_number]
                pa = self.mesh.permeability[global_adj]
                permeability_adj = pa[:, direction_number]
                flow_rate = ((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj)/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()
                effective_permeability.append(center_distance_walls*flow_rate/(area*self.number_faces_coarse_face[direction_number]))

            self.effective_permeability.append(effective_permeability)
            self.effective_porosity.append(effective_porosity)
