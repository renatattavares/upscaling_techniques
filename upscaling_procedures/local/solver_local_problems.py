import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Solver:

    def solve_local_problems(self, transmissibility, source):

        i = self.coarse_volume
        print("Solving local problem {0}".format(i))
        transmissibility = lil_matrix.tocsr(transmissibility)
        source = lil_matrix.tocsr(source)

        pressure = spsolve(transmissibility,source)

        return pressure

        # if np.array_equal(self.direction, self.x) is True:
        #     self.pressure_x = np.append(pressure, self.pressure_x)
        # elif np.array_equal(self.direction, self.y) is True:
        #     self.pressure_y = np.append(pressure, self.pressure_y)
        # elif np.array_equal(self.direction, self.z) is True:
        #     self.pressure_z = np.append(pressure, self.pressure_z)
