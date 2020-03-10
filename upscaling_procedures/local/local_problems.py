"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import copy
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

        # Solve local problems
        self.solve_local_problems()

    def preprocess_mesh(self):

        print('\nPre-processing mesh with IMPRESS...')
        start = time.time()
        self.mesh = impress(self.mesh_file, dim = 3) # IMPRESS' object
        self.coarse = self.mesh.coarse
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

    def set_simulation_variables(self):

        if self.mode is 'auto':
            self.mesh.permeability[:] = np.array([1,1,1])
            self.mesh.porosity[:] = 1

        elif self.mode is 'integrated':
            # self.permeability[:,0] = self.porosity + 5
            # self.permeability[:,1] = self.porosity + 5
            # self.permeability[:,2] = self.porosity + 5
            self.mesh.permeability[:] = self.permeability
            self.mesh.porosity[:] = self.porosity


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

    def solver(self, transmissibility, source):

        i = self.coarse_volume
        #print("Solving local problem {0}".format(i))
        transmissibility = lil_matrix.tocsr(transmissibility)
        source = lil_matrix.tocsr(source)

        pressure = spsolve(transmissibility,source)

        return pressure

    def solve_local_problems(self):

        self.pressure_x = np.array([])
        self.pressure_y = np.array([])
        self.pressure_z = np.array([])

        for i in range(self.number_coarse_volumes):
            print('\nSolving local problem {}'.format(i))
            self.coarse_volume = i
            general_transmissibility = self.assembly_local_problem()
            for j in self.direction_string:
                transmissibility, source = self.set_boundary_conditions(j, general_transmissibility)
                pressure = self.solver(transmissibility, source)

                if j == 'x':
                    self.pressure_x = np.append(self.pressure_x, pressure)
                elif j == 'y':
                    self.pressure_y = np.append(self.pressure_y, pressure)
                elif j == 'z':
                    self.pressure_z = np.append(self.pressure_z, pressure)
