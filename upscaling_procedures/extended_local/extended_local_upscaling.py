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
        #self.upscale_permeability_porosity(0)


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

    def upscale_permeability_porosity(self, cv):

        print('Upscaling of local problem {}'.format(cv))
        volume = self.length_elements[0]*self.length_elements[1]*self.length_elements[2]
        general_transmissibility = self.assembly_extended_local_problem(cv)
        info = []
        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(cv)
        total_volume = volume*len(volumes_global_ids)
        porosity = self.mesh.porosity[volumes_global_ids]
        effective_porosity = (volume*porosity).sum()/total_volume

        for direction in self.direction_string:
            #print('in {} direction'.format(direction))
            self.direction = direction
            transmissibility, source, local_wall = self.set_boundary_conditions(general_transmissibility, cv)
            pressures = self.solver(transmissibility, source)
            direction_number = self.directions_numbers.get(direction)
            center_distance_walls = self.center_distance_walls[cv][direction_number]
            area = self.areas[direction_number]
            global_wall = self.wall(cv)
            local_wall = self.global_to_local_id(volumes_global_ids, global_wall)

            center_wall = self.mesh.volumes.center[global_wall]
            global_adj = self.identify_adjacent_volumes_to_wall(cv, global_wall)
            center_adj = self.mesh.volumes.center[global_adj]
            local_adj = self.global_to_local_id(volumes_global_ids, global_adj)
            pressure_wall = pressures[local_wall]
            pressure_adj = pressures[local_adj]
            pw = self.mesh.permeability[global_wall]
            permeability_wall = pw[:, direction_number]
            pa = self.mesh.permeability[global_adj]
            permeability_adj = pa[:, direction_number]
            flow_rate = ((((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj))*self.areas[direction_number])/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()
            pressure_averaged_gradient = self.pressure_averaged_gradient(pressures, local_wall)
            print(pressure_averaged_gradient)
            info.append(center_distance_walls*flow_rate/(area*len(global_wall)*pressure_averaged_gradient))

        info.append(effective_porosity)

        return info

    def wall(self, coarse_volume):

        fine_volumes = self.coarse.elements[coarse_volume].volumes.father_id[:]
        direction = self.directions_dictionary.get(self.direction)
        direction_number = self.directions_numbers.get(self.direction)
        volumes_centers = self.mesh.volumes.center[fine_volumes]
        coords = volumes_centers[:, direction_number]
        smallest = np.amin(coords)
        wall = np.where(coords == smallest)[0]

        return wall # Global ids

    def identify_adjacent_volumes_to_wall(self, coarse_volume, global_wall):

        fine_volumes = self.coarse.elements[coarse_volume].volumes.father_id[:]
        inter, delete_this, discard = np.intersect1d(fine_volumes, global_wall, return_indices = True)
        fine_volumes = np.delete(fine_volumes, delete_this)
        direction = self.directions_dictionary.get(self.direction)
        direction_number = self.directions_numbers.get(self.direction)
        volumes_centers = self.mesh.volumes.center[fine_volumes]
        coords = volumes_centers[:, direction_number]
        smallest = np.amin(coords)
        adjacent_volumes = np.where(coords == smallest)[0]

        return adjacent_volumes

    def pressure_averaged_gradient(self, pressures, local_wall):

        volume = self.length_elements[0]*self.length_elements[1]*self.length_elements[2]
        fine_volumes = self.coarse.elements[coarse_volume].volumes.father_id[:]
        direction = self.directions_dictionary.get(self.direction)
        direction_number = self.directions_numbers.get(self.direction)
        volumes_centers = self.mesh.volumes.center[fine_volumes]
        coords = volumes_centers[:, direction_number]
        greatest = np.amax(coords)
        end_global_wall = np.where(coords == greatest)[0]
        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(cv)
        end_local_wall = self.global_to_local_id(volumes_global_ids, end_global_wall)
        end_wall_pressures = pressures[end_local_wall]
        wall_pressures = pressures[local_wall]
        volumes = np.arange(len(end_wall_pressures))
        volumes = np.full_like(volumes, volume)
        wall_averaged_pressure = volumes*wall_pressures/volumes.sum()
        end_wall_pressures = volumes*end_wall_pressures/volumes.sum()
        pressure_averaged_gradient = end_wall_pressures - wall_averaged_pressure

        return pressure_averaged_gradient
