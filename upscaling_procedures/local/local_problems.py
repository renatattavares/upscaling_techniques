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

        print('\n##### Treatment of local problems #####')

        # Preprocessing mesh with IMPRESS
        self.mesh_file = mesh_file
        self.preprocess_mesh()

        print('\nMesh informations will be set automatically')

        mesh_info_file = 'mesh_info.yml'
        with open(mesh_info_file, 'r') as file:
            data = yaml.safe_load(file)

        self.number_elements_x_direction = data['x']
        self.number_elements_y_direction = data['y']
        self.number_elements_z_direction = data['z']
        self.number_elements = np.array([self.number_elements_x_direction, self.number_elements_y_direction, self.number_elements_z_direction])
        self.length_elements = np.array([1,1,1])
        self.length_elements_x_direction, self.length_elements_y_direction, self.length_elements_z_direction = 1, 1, 1
        self.mesh.permeability[:] = np.array([1,1,1])
        self.mesh.porosity[:] = 1


        # Read boundary condition chosen
        self.boundary_condition_type = 1

        # Setting variables and informations
        self.set_coordinate_system()
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config())
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

        mesh_info_file = 'mesh_info.yml'
        with open(mesh_info_file, 'r') as file:
            data = yaml.safe_load(file)

        self.direction_string = np.array([])

        if data['x'] != 1:
            self.direction_string = np.append(self.direction_string, 'x')
        if data['y'] != 1:
            self.direction_string = np.append(self.direction_string, 'y')
        if data['z'] != 1:
            self.direction_string = np.append(self.direction_string, 'z')

        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }

    def solver(self, transmissibility, source):

        i = self.coarse_volume
        #print("Solving local problem {0}".format(i))
        transmissibility = lil_matrix.tocsr(transmissibility)
        source = lil_matrix.tocsr(source)

        pressure = spsolve(transmissibility,source)

        return pressure

    def solve_local_problems(self):

        self.pressure_x = []
        self.pressure_y = []
        self.pressure_z = []

        for i in range(self.number_coarse_volumes):
            print('\nSolving local problem {} '.format(i))
            self.coarse_volume = i
            general_transmissibility = self.assembly_local_problem()

            for j in self.direction_string:
                print('{} direction'.format(j))
                transmissibility, source = self.set_boundary_conditions(j, general_transmissibility)
                pressure = self.solver(transmissibility, source)

                if j == 'x':
                    self.pressure_x.append(pressure)
                elif j == 'y':
                    self.pressure_y.append(pressure)
                elif j == 'z':
                    self.pressure_z.append(pressure)

                print('Pressure stored')
