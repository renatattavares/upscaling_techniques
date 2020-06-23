import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class Solver:
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
            general_transmissibility = self.assembly_local_problem()
            p = []
            w = []

            for direction in self.direction_string:
                print('In {} direction'.format(direction))
                self.direction = direction
                transmissibility, source, correct_volumes_group_1 = self.set_boundary_conditions(general_transmissibility)
                p.append(self.solver(transmissibility, source))
                w.append(correct_volumes_group_1)

            self.pressure.append(p)
            self.walls.append(w)

    def solve_permeability_local_upscaling(self):
        pass

    def solve_porosity_local_upscaling(self):
        pass
