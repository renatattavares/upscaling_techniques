"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import numpy as np
import multiprocessing as mp
from upscaling_procedures.local.local_problems import LocalProblems
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ParallelLocalProblems(LocalProblems):

    def __init__(self, mesh_file = None, dataset = None):

        print('\n##### Treatment of local problems #####')

        if mesh_file is None:
            print('\nMesh informations will be accessed from {} dataset'.format(dataset))
            self.mode = 'integrated'
            self.porosity, self.permeability = read_dataset(dataset)
            self.mesh_file = 'generated_mesh.h5m'

        else:
            print('\nMesh informations will be set automatically')
            self.mode = 'auto'
            self.mesh_file = mesh_file

        # Read boundary condition chosen
        self.boundary_condition_type = 1

        # Preprocessing mesh with IMPRESS
        self.preprocess_mesh()

        # Setting variables and informations
        self.set_simulation_variables()
        self.set_coordinate_system()
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config())

        self.pressure_x = np.array([])
        self.pressure_y = np.array([])
        self.pressure_z = np.array([])

        self.distribute_data()
        self.create_processes()

        # Solve local problems
        #self.solve_local_problems()
    def solve_local_problems(self, coarse_volume):

        p = []

        print('\nSolving local problem {}'.format(coarse_volume))
        self.coarse_volume = coarse_volume
        general_transmissibility = self.assembly_local_problem()
        for j in self.direction_string:
            transmissibility, source = self.set_boundary_conditions(j, general_transmissibility)
            pressure = self.solver(transmissibility, source)
            p.append(pressure)

        return p

    def solve_local_problems_parallel(self, coarse_volumes, queue):

        for cv in coarse_volumes:
            p = self.solve_local_problems(cv)

        queue.put(p)

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

    def create_processes(self):

        pressures = []
        process_number = len(self.distribution)

        queues = [mp.Queue() for i in range(process_number)]
        processes = [mp.Process(target = self.solve_local_problems_parallel, args = (i,q,)) for i, q in zip(self.distribution, queues)]

        for p in processes:
            p.start()

        for p, q in zip(processes, queues):
            pressures.append(q.get())
            p.join()

        self.pressures = pressures
