"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
#import xlsxwriter
from upscaling_procedures.local.solver_local_problems import Solver
from upscaling_procedures.local.boundary_conditions import BoundaryConditions
from upscaling_procedures.local.assembly import Assembly
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

class LocalProblems(BoundaryConditions, Solver, Assembly):

    def __init__(self, mesh_file = None, boundary_condition_type = None):

        print('\n##### Treatment of local problems #####')

        # Preprocessing mesh with IMPRESS
        print('\nPre-processing mesh with IMPRESS...')
        start = time.time()
        self.mesh = impress(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Setting variables
        print('\nAccessing coarsening informations from IMPRESS and setting important variables...')
        self.mesh.permeability[:] = np.array([1,1,1])
        self.mesh.porosity[:] = 1
        self.mesh.equivalent_permeability[:] = 0
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.direction_string = np.array(['x', 'y', 'z'])
        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config())
        self.coarse = self.mesh.coarse
        self.boundary_condition_type = boundary_condition_type
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.pressure_gradient = 1
        self.pressure_x = np.array([])
        self.pressure_y = np.array([])
        self.pressure_z = np.array([])
        self.correct_volumes_x = np.array([], dtype = int)
        self.correct_volumes_y = np.array([], dtype = int)
        self.correct_volumes_z = np.array([], dtype = int)

        self.run()

    def run(self):
        for i in self.direction_string:
            print('\n##### Local problems in {} direction #####'.format(i))

            for j in range(self.number_coarse_volumes):
                self.coarse_volume = j
                transmissibility = self.assembly_local_problem()
                transmissibility, source = self.set_boundary_conditions(i,transmissibility)
                pressure = self.solve_local_problems(transmissibility, source)
                print('\n')

                if i == 'x':
                    self.pressure_x = np.append(self.pressure_x, pressure)
                elif i == 'y':
                    self.pressure_y = np.append(self.pressure_y, pressure)
                elif i == 'z':
                    self.pressure_z = np.append(self.pressure_z, pressure)
