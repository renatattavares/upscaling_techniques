# LOCAL UPSCALLING OF STRUCTURED MESHES IN HOMOGENEOS MEDIA
import sys
sys.path.append('/elliptic_case')
from static.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as preprocessor
import numpy as np
import time
import pdb
from math import pi
from pymoab import rng
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
#import xlsxwriter

M = preprocessor('mesh/25.h5m', dim = 3)
dx, dy, dz = 1, 1, 1
nx, ny, nz = np.cbrt(len(M.volumes)), np.cbrt(len(M.volumes)), np.cbrt(len(M.volumes))
cx, cy, cz = 5, 5, 5 # Coarsening ratio
rx, ry, rz = (nx/cx), (ny/cy), (nz/cz)
num_elements = len(M.volumes)
num_elements_coarse = rx*ry*rz
num_elements_coarse = np.int64(num_elements_coarse)

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return ((c1-c2)**2).sum()

print("Setting the permeability tensor")
M.permeability[:] = np.array([1, 1, 1])

for i in range(len(M.coarse.elements)):
    print("Assembly of coarse volume {0}".format(i))
    i = np.int64(i)
    start = time.time()
    adj = M.coarse.elements[i].volumes.bridge_adjacencies(M.coarse.elements[i].volumes.all, 2, 3) # IDs locais
    perm = M.permeability[np.int64(M.coarse.elements[i].volumes.global_id[np.int64(M.coarse.elements[i].volumes.all)])]
    center = M.coarse.elements[i].volumes.center[np.int64(M.coarse.elements[i].volumes.all)]
    coef = lil_matrix(num_elements_coarse, num_elements_coarse)
    for b in range(num_elements_coarse):
        adjacencies = adj[b] # Array de IDs locais
        for c in range(len(adjacencies)):
            id = np.array( [adjacencies[c]],  dtype= np.int64)
            coef[b,id[0]] = equiv_perm(perm[b], perm[id])/1
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

    print("Solving the problem of coarse volume {0}".format(i))
    start = time.time()
    coef = lil_matrix.tocsr(coef)
    q = lil_matrix.tocsr(q)
    P_coarse_volume = spsolve(coef,q)
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Storing results of coarse volume {0}".format(i))
    start = time.time()
    M.coarse.elements[i].pressure_coarse[:] = P_coarse_volume
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    total_flow = 0.0
    flow_rate = 0.0
    for v in range(rx*ry):
        flow_rate =  + equiv_perm(perm[v], perm[v+rx*ry])*area*(M.coarse.elements[i].pressure_coarse[v]-M.coarse.elements[i].pressure_coarse[v+rx*ry])
        total_flow = total_flow + flow_rate

    permeability_coarse = total_flow/((area*rx*ry)*(M.coarse.elements[i].pressure_coarse[v]-M.coarse.elements[i].pressure_coarse[v+rx*ry]))
    print(permeability_coarse)

print("Assembly of upscaling")
start = time.time()
coef = lil_matrix((len(M.coarse.elements), len(M.coarse.elements)), dtype=np.float_)

for i in range(len(M.coarse.elements)):
    #M.coarse_volumes[i].permeability_coarse[:] = permeability_coarse
    #perm = M.coarse_volumes[i].permeability_coarse[:]
    adj = M.coarse.elements[i].faces.coarse_neighbors
    for j in range(len(adj)):
        id = np.array(adj[j],  dtype= np.int)
        coef[i,id] = equiv_perm(1, 1)/25
    coef[i,i] = (-1)*coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting boundary conditions of coarse mesh")
start = time.time()
q = lil_matrix((len(M.coarse.elements), 1), dtype=np.float_)
coef[0:25] = 0
q [0:25] = 500
coef[100:125] = 0
for r in range(25):
    coef[r,r] = 1
    coef[r+100,r+100] = 1
end = time.time()
print("This step lasted {0}s".format(end-start))

'''
workbook = xlsxwriter.Workbook('correct.xlsx')
worksheet = workbook.add_worksheet()
matrix = lil_matrix.toarray(coef)

row = 0
col = 0

for row in range(125):
  for col in range(125):
    worksheet.write(row, col, matrix[row][col])

workbook.close()
'''

print("Solving the problem")
start = time.time()
coef = lil_matrix.tocsr(coef)
q = lil_matrix.tocsr(q)
P = spsolve(coef,q)
end = time.time()
print("This step lasted {0}s".format(end-start))


for cvolume,index in zip(M.coarse.elements,range(len(P))):
    M.pressure[cvolume.volumes.global_id[cvolume.volumes.all]] = P[index]

print("Printing results")
M.core.print()
