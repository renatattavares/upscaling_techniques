import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Solver:

    def solve_local_problems(self):

        for i in range(self.number_coarse_volumes):
            print("Solving local problem of coarse volume {0}".format(i))
            transmissibility = lil_matrix.tocsr(self.transmissibilities[i])
            source = lil_matrix.tocsr(self.sources[i])

            if np.array_equal(self.direction, self.x) is True:
                self.pressure_x = np.array([])
                self.pressure_x = np.append(self.pressure_x, spsolve(transmissibility,source))
            elif np.array_equal(self.direction, self.y) is True:
                self.pressure_y = np.array((len(self.mesh.volumes), 1))
                self.pressure_y = np.append(self.pressure_y, spsolve(transmissibility,source))
            elif np.array_equal(self.direction, self.z) is True:
                self.pressure_z = np.array((len(self.mesh.volumes), 1))
                self.pressure_z = np.append(self.pressure_z, spsolve(transmissibility,source))
