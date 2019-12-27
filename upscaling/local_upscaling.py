# LOCAL UPSCALLING OF STRUCTURED MESHES IN HOMOGENEOS MEDIA
import numpy as np
import time
import pdb
import xlsxwriter
from math import pi
from pymoab import rng
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
from preprocessor import M

dx, dy, dz = 1, 1, 1
nx, ny, nz = 25, 25,25
cx, cy, cz = 5, 5, 5
rx, ry, rz = 5, 5, 5
num_elements = nx*ny*nz
num_elements_coarse = rx*ry*rz

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return ((c1-c2)**2).sum()

print("Setting the permeability tensor")
M.permeability[:] = 1 #np.array([[1, 0],[2, 0],[3, 0]])
area = dx*dy

for i in range(len(M.coarse_volumes)):
    print("Assembly of coarse volume {0}".format(i))
    start = time.time()
    adj = M.coarse_volumes[i].volumes.bridge_adjacencies(M.coarse_volumes[i].volumes.all, 2, 3) # IDs locais
    perm = M.permeability[M.coarse_volumes[i].volumes.global_id[M.coarse_volumes[i].volumes.all]]
    center = M.coarse_volumes[i].volumes.center[M.coarse_volumes[i].volumes.all]
    coef = lil_matrix((num_elements_coarse, num_elements_coarse), dtype=np.float_)
    for b in range(num_elements_coarse):
        adjacencies = adj[b] # Array de IDs locais
        for c in range(len(adjacencies)):
            id = np.array( [adjacencies[c]],  dtype= np.int)
            coef[b,id] = 1
        coef[b,b] = (-1)*coef[b].sum()
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Setting boundary conditions of coarse volume {0}".format(i))
    start = time.time()
    q = lil_matrix((num_elements_coarse, 1), dtype=np.float_)
    coef[0:rx*ry] = 0
    q [0:rx*ry] = 500
    coef[(num_elements_coarse)-(rx*ry):num_elements_coarse] = 0
    for r in range(rx*ry):
        coef[r,r] = 1
        coef[r+(num_elements_coarse)-(rx*ry),r+(num_elements_coarse)-(rx*ry)] = 1
    end = time.time()
    print("This step lasted {0}s".format(end-start))
