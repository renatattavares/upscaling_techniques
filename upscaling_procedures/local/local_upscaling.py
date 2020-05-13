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
        self.upscale_permeability()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))


    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        area = 1
        for i in range(self.number_coarse_volumes):
            for j in self.direction_string:
                direction = j

                if direction == 'x':
                    local_wall = self.wall_x[i]
                    pressures = self.pressure_x[i]
                    center_distance_walls = self.center_distance_walls_x[i]
                elif direction == 'y':
                    local_wall = self.wall_y[i]
                    pressures = self.pressure_y[i]
                    center_distance_walls = self.center_distance_walls_y[i]
                elif direction == 'z':
                    local_wall = self.wall_z[i]
                    pressures = self.pressure_z[i]
                    center_distance_walls = self.center_distance_walls_z[i]

                global_wall = self.coarse.elements[i].volumes.global_id[local_wall]
                local_adj = self.identify_adjacent_volumes_to_wall(i, local_wall)
                global_adj = self.coarse.elements[i].volumes.global_id[local_adj]
                center_wall = self.coarse.elements[i].volumes.center[local_wall]
                center_adj = self.coarse.elements[i].volumes.center[local_adj]
                pressure_wall = pressures[local_wall]
                pressure_adj = pressures[local_adj]
                permeability_wall = self.get_absolute_permeabilities(direction, global_wall)
                permeability_adj = self.get_absolute_permeabilities(direction, global_adj)

                flow_rate = ((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj)/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()

                effective_permeability = center_distance_walls*flow_rate/(area*self.number_faces_coarse_face)

                print(effective_permeability)

    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.global_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.ones(len(porosity))
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()

    def get_absolute_permeabilities(self, direction, global_ids):

        if direction == 'x':
            permeability = self.mesh.permeability[global_ids]
            permeability = permeability[:,0]
        elif direction == 'y':
            permeability = self.mesh.permeability[global_ids]
            permeability = permeability[:,1]
        elif direction == 'z':
            permeability = self.mesh.permeability[global_ids]
            permeability = permeability[:,2]

        return permeability
