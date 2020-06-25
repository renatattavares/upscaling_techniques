"""
Module of local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import yaml
import numpy as np
import multiprocessing as mp
from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.visualize import Visualize
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ParallelLocalUpscaling(LocalUpscaling, Visualize):

    def __init__(self, mesh_file = None, dataset = None):

        initial_time = time.time()

        super().__init__(mesh_file, dataset)
        # Upscale in parallel
        self.distribute_data()
        self.create_processes()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def distribute_data(self):

        core_number = mp.cpu_count() # Cores in the computer
        number_problems_in_process = int(self.number_coarse_volumes/core_number) # Average number of local problems to be solved by one process
        rest = self.number_coarse_volumes % core_number
        distribution = []
        coarse_volumes = np.arange(self.number_coarse_volumes)

        for i in range(core_number):
            if i is core_number-1:
                distribution.append(coarse_volumes[number_problems_in_process*i:number_problems_in_process*(i+1)+rest])
            else:
                distribution.append(coarse_volumes[number_problems_in_process*i:number_problems_in_process*(i+1)])

        self.distribution = distribution

    def upscale_permeability_porosity(self, cv):

        print('Upscaling of local problem {}'.format(cv))
        self.coarse_volume = cv
        volume = self.length_elements[0]*self.length_elements[1]*self.length_elements[2]
        general_transmissibility = self.assembly_local_problem()
        info = []
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
            flow_rate = ((((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj))*self.areas[direction_number])/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()
            info.append(center_distance_walls*flow_rate/(area*self.number_faces_coarse_face[direction_number]))

        info.append(effective_porosity)

        return info

    def upscale_permeability_porosity_parallel(self, coarse_volumes, queue):

        info = []

        for cv in coarse_volumes:
            info.append(self.upscale_permeability_porosity(cv))

        queue.put(info)

    def create_processes(self):

        self.info = []
        process_number = len(self.distribution)

        queues = [mp.Queue() for i in range(process_number)]
        processes = [mp.Process(target = self.upscale_permeability_porosity_parallel, args = (i,q,)) for i, q in zip(self.distribution, queues)]

        for p in processes:
            p.start()

        for p, q in zip(processes, queues):
            self.info.append(q.get())
            p.join()
