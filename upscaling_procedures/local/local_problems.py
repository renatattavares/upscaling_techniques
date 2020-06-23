"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import copy
import yaml
#import xlsxwriter
import numpy as np
from upscaling_procedures.local.solver import Solver
from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.assembly import Assembly
from upscaling_procedures.local.mesh_geometry import MeshGeometry
from upscaling_procedures.local.boundary_conditions import BoundaryConditions
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class LocalProblems(MeshGeometry, Assembly, BoundaryConditions, Solver):

    def __init__(self, mesh_file = None, dataset = None):

        print('\n##### Treatment of local problems #####')

        if mesh_file is None:
            self.mode = 'integrated'

        else:
            self.mode = 'auto'
            self.mesh_file = mesh_file

        # Preprocessing mesh with IMPRESS
        self.preprocess_mesh()

        # Setting variables and informations
        self.boundary_condition_type = 1 # Fixed pressure
        self.set_coordinate_system()
        self.set_simulation_variables()

        # Preparing to solve local problems
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config())

        # Solve local problems
        self.solve_local_problems()

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

        if self.mode is 'auto':
            self.mesh.permeability[:] = np.array([1,1,1])
            self.mesh.porosity[:] = 1

            mesh_info_file = 'mesh_info.yml'

            with open(mesh_info_file, 'r') as file:
                data = yaml.safe_load(file)

            self.number_elements_x_direction = data['Elements']['x']
            self.number_elements_y_direction = data['Elements']['y']
            self.number_elements_z_direction = data['Elements']['z']

            self.length_elements_x_direction = data['Lenght']['x']
            self.length_elements_y_direction = data['Lenght']['y']
            self.length_elements_z_direction = data['Lenght']['z']

        elif self.mode is 'integrated':
            self.mesh_file = 'mesh/dataset_mesh.h5m'
            self.porosity, self.permeability, self.number_elements, self.length_elements = read_dataset(dataset)

            self.number_elements_x_direction = self.number_elements[0]
            self.number_elements_y_direction = self.number_elements[1]
            self.number_elements_z_direction = self.number_elements[2]

            self.length_elements_x_direction = self.length_elements[0]
            self.length_elements_y_direction = self.length_elements[1]
            self.length_elements_z_direction = self.length_elements[2]
