import numpy as np
import multiprocessing as mp
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling

class FluxTest(ParallelLocalUpscaling):

    def upscale_permeability(self, coarse_volume, pressure):

        area = 1
        i = coarse_volume
        effective_permeability = []
        flux = []
        info = []

        for j in self.direction_string:
            direction = j
            if direction == 'x':
                local_wall = self.get_wall(i, 0)
                pressures = pressure[0]
                center_distance_walls = self.center_distance_walls_x[i]
            elif direction == 'y':
                local_wall = self.get_wall(i, 1)
                pressures = pressure[1]
                center_distance_walls = self.center_distance_walls_y[i]
            elif direction == 'z':
                local_wall = self.get_wall(i, 2)
                pressures = pressure[2]
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

            flux.append(flow_rate)
            effective_permeability.append(center_distance_walls*flow_rate/(area*self.number_faces_coarse_face))

        info.append(effective_permeability)
        info.append(flux)

        return info

    def upscale_permeability_parallel(self, coarse_volumes, queue):

        info = []

        for cv in coarse_volumes:
            p = self.solve_local_problems(cv)
            info.append(self.upscale_permeability(cv, p))
            print('Done with coarse volume {}'.format(cv))

        queue.put(info)

    def create_processes(self):

        info = []
        process_number = len(self.distribution)

        queues = [mp.Queue() for i in range(process_number)]
        processes = [mp.Process(target = self.upscale_permeability_parallel, args = (i,q,)) for i, q in zip(self.distribution, queues)]

        for p in processes:
            p.start()

        for p, q in zip(processes, queues):
            info.append(q.get())
            p.join()

        self.info = info
