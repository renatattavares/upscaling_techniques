"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import copy
import yaml
#import xlsxwriter
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.assembly import Assembly
from upscaling_procedures.local.boundary_conditions import BoundaryConditions
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class LocalProblems(Assembly, BoundaryConditions):

    def __init__(self, mesh_file = None, dataset = None):

        print('\n##### LocalProblems class initialized #####')

        if mesh_file is None:
            self.mesh_file = 'mesh/dataset_mesh.h5m'
            porosity, permeability, self.number_elements, self.length_elements = read_dataset(dataset, self.mesh_file)
            # Preprocessing mesh with IMPRESS
            self.preprocess_mesh()
            self.mesh.porosity[:] = porosity
            self.mesh.permeability[:] = permeability

        else:
            self.mode = 'auto'
            self.mesh_file = mesh_file
            # Preprocessing mesh with IMPRESS
            self.preprocess_mesh()
            self.set_simulation_variables()


        # Setting variables and informations
        self.boundary_condition_type = 1 # Fixed pressure
        self.set_coordinate_system()

        # Preparing to solve local problems
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config())

        # Solve local problems
        #self.solve_local_problems()

    def preprocess_mesh(self):

        print('\nPre-processing mesh with IMPRESS...')
        start = time.time()
        self.mesh = impress(self.mesh_file, dim = 3) # IMPRESS' object
        self.coarse = self.mesh.coarse
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

    def set_coordinate_system(self):
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.direction_string = np.array(['x', 'y', 'z'])
        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }
        self.directions_numbers = {
            'x': 0,
            'y': 1,
            'z': 2
            }

    def set_simulation_variables(self):

        self.mesh.permeability[:] = np.array([1,1,1])
        self.mesh.porosity[:] = 1

        mesh_info_file = 'mesh_info.yml'

        with open(mesh_info_file, 'r') as file:
            data = yaml.safe_load(file)

        number_elements = data['Elements']
        lenght_elements = data['Length']
        self.number_elements = np.array([], dtype = int)
        self.length_elements = np.array([], dtype = int)

        for number, lenght in zip(number_elements.values(), lenght_elements.values()):
            self.number_elements = np.append(self.number_elements, number)
            self.length_elements = np.append(self.length_elements, lenght)

    def solver(self, transmissibility, source):

        transmissibility = lil_matrix.tocsr(transmissibility)
        source = lil_matrix.tocsr(source)
        pressure = spsolve(transmissibility,source)

        return pressure

    def solve_local_problems(self):
        """
        Implementation of solver for LocalProblem class
        """
        self.pressure = []
        self.walls = []

        for coarse_volume in range(self.number_coarse_volumes):
            print('\nSolving local problem {}'.format(coarse_volume))
            self.coarse_volume = coarse_volume
            general_transmissibility = self.assembly_local_problem(coarse_volume)
            p = []
            w = []

            for direction in self.direction_string:
                print('In {} direction'.format(direction))
                self.direction = direction
                transmissibility, source, correct_volumes_group_1 = self.set_boundary_conditions(general_transmissibility, coarse_volume)
                self.transmissibility = transmissibility
                self.source = source
                p.append(self.solver(transmissibility, source))
                w.append(correct_volumes_group_1)

            self.pressure.append(p)
            self.walls.append(w)
