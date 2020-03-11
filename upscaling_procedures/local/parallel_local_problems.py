"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import numpy as np
from multiprocessing import Process, Queue
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

        self.create_processes()

        # Solve local problems
        #self.solve_local_problems()
    def solve_local_problems(self, coarse_volume, queue):

        print('\nSolving local problem {}'.format(coarse_volume))
        self.coarse_volume = coarse_volume
        general_transmissibility = self.assembly_local_problem()
        for j in self.direction_string:
            transmissibility, source = self.set_boundary_conditions(j, general_transmissibility)
            pressure = self.solver(transmissibility, source)

        #queue.put(pressure)

    def distribute_data(self, q):
        pass

    def create_processes(self):

        process_number = len(self.mesh.coarse.elements)

        #Processes' arguments
        queues = [Queue() for i in range(process_number)]
        coarse_volumes = np.arange(process_number)

        # Create and initialize processes
        processes = [Process(target = self.solve_local_problems, args = (i, q,)) for i, q in zip(coarse_volumes, queues)]

        for p in processes:
            p.start()

        for p in processes:
            p.join()
