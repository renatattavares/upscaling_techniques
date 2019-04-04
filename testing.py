# LOCAL UPSCALLING OF STRUCTURED MESHES IN HOMOGENEOUS MEDIA

import numpy as np
import time
import pdb
import xlsxwriter
from math import pi
from pymoab import rng, types
from boundary_conditions import BoundaryConditions
from flow_based_permeability import FlowBasedPermeability
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
from preprocessor import M

dx, dy, dz = 1, 1, 1
nx, ny, nz = 5, 5,5
cx, cy, cz = 1, 1, 1
rx, ry, rz = 5, 5, 5
num_elements = nx*ny*nz
num_elements_coarse = rx*ry*rz

def equiv_perm(k1, k2):
    return (2*k1*k2)/(k1 + k2)

def centroid_dist(c1, c2):
    return ((c1-c2)**2).sum()

print("Setting the permeability")
#M.permeability[:] = 1
permeability = np.array([[1,0,0],[0,100,0],[0,0,200]])
flux_direction = ['x', 'z']
area = dx*dy

# for d in flux_direction:
#     print("Upscaling in {0} direction".format(d))

for i in range(len(M.coarse_volumes)):
    print("Assembly of coarse volume {0}".format(i))
    start = time.time()
    adj = M.coarse_volumes[i].volumes.bridge_adjacencies(M.coarse_volumes[i].volumes.all, 2, 3) # IDs locais
    #perm = M.permeability[M.coarse_volumes[i].volumes.global_id[M.coarse_volumes[i].volumes.all]]
    center = M.coarse_volumes[i].volumes.center[M.coarse_volumes[i].volumes.all]
    coef = lil_matrix((num_elements_coarse, num_elements_coarse), dtype=np.float_)
    for b in range(num_elements_coarse):
        adjacencies = adj[b]
        for c in range(len(adjacencies)):
            id = np.array([adjacencies[c]],  dtype= np.int)
            direction = (center[b]-center[id])**2
            perm = np.linalg.norm(np.dot(direction, permeability))
            coef[b,id] = perm/centroid_dist(center[b], center[id])
            print(coef[b,id])
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
    M.coarse_volumes[i].pressure_coarse[:] = P_coarse_volume
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Calculating effective permeability")
    start = time.time()

    total_flow = 0.0
    flow_rate = 0.0

    for v in range(rx*ry):
        adjacencies = adj[v]
        for w in range(len(adjacencies)):
            id = np.array([adjacencies[w]],  dtype= np.int)
            direction = (center[v]-center[id])**2
            perm = np.linalg.norm(np.dot(direction, permeability))
            flow_rate =  + perm*area*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[id])
            total_flow = total_flow + flow_rate

    permeability_coarse = total_flow/((area*rx*ry)*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[id]))
    print(permeability_coarse)

    end = time.time()
    print("This step lasted {0}s".format(end-start))


print("Assembly of upscaling")
start = time.time()
coarse_coef = lil_matrix((len(M.coarse_volumes), len(M.coarse_volumes)), dtype=np.float_)

for i in range(len(M.coarse_volumes)):
    #M.coarse_volumes[i].permeability_coarse[:] = permeability_coarse
    #perm = M.coarse_volumes[i].permeability_coarse[:]
    adj = M.coarse_volumes[i].faces.coarse_neighbors
    for j in range(len(adj)):
        id = np.array(adj[j],  dtype= np.int)
        coarse_coef[i,id] = 1/2.5
    coarse_coef[i,i] = (-1)*coarse_coef[i].sum()
end = time.time()
print("This step lasted {0}s".format(end-start))

print("Setting boundary conditions of coarse mesh")
start = time.time()
coarse_coef, coarse_q = bc.coarse_bc(coarse_coef, len(M.coarse_volumes))
end = time.time()
print("This step lasted {0}s".format(end-start))

'''
workbook = xlsxwriter.Workbook('verify.xlsx')
worksheet = workbook.add_worksheet()
matrix = lil_matrix.toarray(coarse_coef)

row = 0
col = 0

for row in range(125):
  for col in range(125):
    worksheet.write(row, col, matrix[row][col])

workbook.close()
'''

print("Solving the problem")
start = time.time()
coarse_coef = lil_matrix.tocsr(coarse_coef)
coarse_q = lil_matrix.tocsr(coarse_q)
coarse_p = spsolve(coarse_coef,coarse_q)
end = time.time()
print("This step lasted {0}s".format(end-start))

for cvolume,index in zip(M.coarse_volumes,range(len(coarse_p))):
    M.pressure[cvolume.volumes.global_id[cvolume.volumes.all]] = coarse_p[index]

print("Printing results")
M.core.print()
