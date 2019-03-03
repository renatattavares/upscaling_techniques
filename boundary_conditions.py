import numpy as np
from scipy.sparse import csr_matrix, lil_matrix

class BoundaryConditions():
    def __init__(self, num_elements_coarse, rx, ry, coef, flux_direction, coarse_mesh):
        self.coef = coef
        self.num_elements_coarse = num_elements_coarse
        self.rx = rx
        self.ry = ry
        self.q = lil_matrix((self.num_elements_coarse, 1), dtype=np.float_)
        self.flux_direction = flux_direction
        self.coarse_mesh = coarse_mesh

        if self.coarse_mesh != 1:
            if self.flux_direction == 'x':
                self.coef, self.q = self.x()
            elif self.flux_direction == 'y':
                pass
            elif self.flux_direction == 'z':
                self.coef, self.q = self.z()
            return self.coef, self.q
        # elif
        #     if self.flux_direction == 'x':
        #         pass
        #         #self.coef, self.q = self.coarse_x()
        #     elif self.flux_direction == 'y':
        #         pass
        #     elif self.flux_direction == 'z':
        #         pass
        #         #self.coef, self.q = self.coarse_z()

    def x(self):
        elements = np.array((), dtype=int)
        elements2 = np.array((), dtype=int)

        for i in range(self.rx*self.ry):
            elements = np.append(elements, ((i+1)-1)*5)
            #print(((i+1)-1)*5)
        for i in range(self.rx*self.ry):
            elements2 = np.append(elements2, ((self.rx-1) + ((i+1)-1)*5))
            #print(4 + ((i+1)-1)*5)

        self.coef[elements] = 0
        self.q [elements] = 500
        self.coef[elements2] = 0
        for r in range(self.rx*self.ry):
            self.coef[elements[r],elements[r]] = 1
            self.coef[elements2[r],elements2[r]] = 1
        return self.coef, self.q

    def y(self):
        pass

    def z(self):
        self.coef[0:rx*ry] = 0
        self.q [0:rx*ry] = 500
        self.coef[(self.num_elements_coarse)-(self.rx*self.ry):self.num_elements_coarse] = 0
        for r in range(self.rx*self.ry):
            coef[r,r] = 1
            coef[r+(self.num_elements_coarse)-(self.rx*self.ry),r+(self.num_elements_coarse)-(self.rx*self.ry)] = 1
        return self.coef, self.q
