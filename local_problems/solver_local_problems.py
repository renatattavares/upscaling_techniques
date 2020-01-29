import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Solver:

    def solve_local_problems(self):

        for i in range(self.number_coarse_volumes):
            print("Solving local problem of coarse volume {0}".format(i))
            transmissibility = lil_matrix.tocsr(self.coarse.elements[i].transmissibility)
            source = lil_matrix.tocsr(self.coarse.elements[i].source)
            global_id_volumes = self.coarse.elements[i].volumes.father_id[:]

            if np.array_equal(self.direction, self.x) == True:
                self.mesh.pressure_x[global_id_volumes] = spsolve(transmissibility,source)
            elif np.array_equal(self.direction, self.y) == True:
                self.mesh.pressure_y[global_id_volumes] = spsolve(transmissibility,source)
            elif np.array_equal(self.direction, self.z) == True:
                self.mesh.pressure_z[global_id_volumes] = spsolve(transmissibility,source)
