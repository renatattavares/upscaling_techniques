import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Solver:

    def solve_local_problems(self):

        if np.array_equal(self.direction, self.x) is True:
            self.pressure_x = np.array([])
        elif np.array_equal(self.direction, self.y) is True:
            self.pressure_y = np.array([])
        elif np.array_equal(self.direction, self.z) is True:
            self.pressure_z = np.array([])

        for i in range(self.number_coarse_volumes):
            print("Solving local problem of coarse volume {0}".format(i))
            transmissibility = lil_matrix.tocsr(self.transmissibilities[i])
            source = lil_matrix.tocsr(self.sources[i])

            pressure = spsolve(transmissibility,source)

            if np.array_equal(self.direction, self.x) is True:
                self.pressure_x = np.append(pressure, self.pressure_x)
            elif np.array_equal(self.direction, self.y) is True:
                self.pressure_y = np.append(pressure, self.pressure_y)
            elif np.array_equal(self.direction, self.z) is True:
                self.pressure_z = np.append(pressure, self.pressure_z)
