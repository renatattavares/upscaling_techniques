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

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))


    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        area = 1
        for i in range(self.number_coarse_volumes):
            for j in self.direction_string:
                direction = self.directions_dictionary.get(j)

                if direction is 'x':
                    wall = self.wall_x[i]
                    pressures = self.pressure_x
                elif direction is 'y':
                    wall = self.wall_y[i]
                    pressures = self.pressure_y
                elif direction is 'z':
                    wall = self.wall_z[i]
                    pressures = self.pressure_z

                global_wall = self.coarse.elements[i].volumes.global_id[wall]
                global_adj = identify_adjacent_volumes_to_wall(wall)
                center_wall = self.coarse.elements[i].volumes.center[wall]
                center_adj = self.coarse.elements[i].volumes.center[adjacent_volumes]
                pressure_wall = pressures[global_wall]
                pressure_adj = pressure[global_adj]
                permeability_wall = get_absolute_permeabilities(direction, global_wall)
                permeability_adj = get_absolute_permeabilities(direction, global_adj)

                flow_rate = (2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj)/np.linalg.norm(center_wall - center_adj, axis = 1).sum()

                effective_permeability = 4*self.flow_rate/(area*self.number_faces_coarse_face)

    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.global_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.ones(len(porosity))
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()

    def identify_adjacent_volumes_to_wall(self, wall):

        adjacent_volumes = self.coarse.elements[i].volumes.bridge_adjacencies(wall, 2, 3).flatten()
        repeated_volumes = np.isin(adjacent, wall, invert = True)
        adjacent_volumes = np.unique(adjacent_volumes[repeated_volumes])
        global_id_adjacent_volumes = self.coarse.elements[i].volumes.global_id[adjacent_volumes]

        return global_id_adjacent_volumes

    def get_absolute_permeabilities(self, direction, global_ids):

        if direction is 'x':
            permeability = self.mesh.permeability[global_ids][0]
        elif direction is 'y':
            permeability = self.mesh.permeability[global_ids][1]
        elif direction is 'z':
            permeability = self.mesh.permeability[global_ids][2]

        return permeability
