import numpy as np
from scipy.sparse import csr_matrix, lil_matrix

class BoundaryConditions():
    def __init__(self, num_elements_coarse, rx, ry, coef, flux_direction):
        self.coef = coef
        self.num_elements_coarse = num_elements_coarse
        self.rx = rx
        self.ry = ry
        self.q = lil_matrix((self.num_elements_coarse, 1), dtype=np.float_)
        self.flux_direction = flux_direction
        self.elements = np.array((), dtype=int)
        self.elements2 = np.array((), dtype=int)

        if self.flux_direction == 'x':
            self.coef, self.q = self.x()
        elif self.flux_direction == 'y':
            pass
        elif self.flux_direction == 'z':
            self.coef, self.q = self.z()
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
        print('Flux direction x')

        for i in range(self.rx*self.ry):
            self.elements = np.append(self.elements, ((i+1)-1)*5)
            #print(((i+1)-1)*5)
        for i in range(self.rx*self.ry):
            self.elements2 = np.append(self.elements2, ((self.rx-1) + ((i+1)-1)*5))
            #print(4 + ((i+1)-1)*5)

        self.coef[self.elements] = 0
        self.q [self.elements] = 500
        self.coef[self.elements2] = 0
        for r in range(self.rx*self.ry):
            self.coef[self.elements[r],self.elements[r]] = 1
            self.coef[self.elements2[r],self.elements2[r]] = 1
        return self.coef, self.q

    def y(self):
        pass

    def z(self):

        self.elements = np.arange(self.rx*self.ry)
        self.elements2 = np.arange(((self.num_elements_coarse)-(self.rx*self.ry)),self.num_elements_coarse)

        self.coef[self.elements] = 0
        self.q [self.elements] = 500
        self.coef[self.elements2] = 0
        for r in range(self.rx*self.ry):
            self.coef[self.elements[r],self.elements[r]] = 1
            self.coef[self.elements2[r],self.elements2[r]] = 1
        return self.coef, self.q

    def coarse_bc(self, coarse_coef, coarses):

        self.coarse_coef = coarse_coef
        self.coarse_q = lil_matrix((coarses, 1), dtype=np.float_)

        for r in range(self.rx*self.ry):
            self.coarse_coef[self.elements[r]] = 0
            self.coarse_q [self.elements[r]] = 500
            self.coarse_coef[self.elements2[r]] = 0
            self.coarse_coef[self.elements[r],self.elements[r]] = 1
            self.coarse_coef[self.elements2[r],self.elements2[r]] = 1

        return self.coarse_coef, self.coarse_q
