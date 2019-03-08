# LOCAL UPSCALLING OF STRUCTURED MESHES IN HOMOGENEOUS MEDIA

import numpy as np
import time
import pdb
import xlsxwriter
from math import pi
from pymoab import rng, types
from boundary_conditions import BoundaryConditions
from scipy.sparse import csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
import mspreprocessor.geoUtil.geoTools as gtool
from mspreprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh

start = time.time()
M = msh("25.h5m", dim = 3)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))

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

print("Setting the permeability")
M.permeability[:] = 1
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
            coef[b,id] = equiv_perm(perm[b], perm[id])/1
        coef[b,b] = (-1)*coef[b].sum()
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Setting boundary conditions of coarse volume {0}".format(i))
    start = time.time()
    flux_direction = 'z'
    bc = BoundaryConditions(num_elements_coarse, rx,ry, coef, flux_direction)
    end = time.time()
    print("This step lasted {0}s".format(end-start))

    print("Solving the problem of coarse volume {0}".format(i))
    start = time.time()
    coef = lil_matrix.tocsr(bc.coef)
    q = lil_matrix.tocsr(bc.q)
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

    for v in range(len(bc.elements)):
        flow_rate =  + equiv_perm(perm[v], perm[v+rx*ry])*area*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+rx*ry])
        total_flow = total_flow + flow_rate

    permeability_coarse = total_flow/((area*rx*ry)*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+rx*ry]))
    #     flow_rate =  + equiv_perm(perm[v], perm[v+1])*area*(M.coarse_volumes[i].pressure_coarse[np.array(v)]-M.coarse_volumes[i].pressure_coarse[np.array(v+1)])
    #     total_flow = total_flow + flow_rate
    #
    # permeability_coarse = total_flow/((area*rx*ry)*(M.coarse_volumes[i].pressure_coarse[0]-M.coarse_volumes[i].pressure_coarse[1]))
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
        local_id = np.array(adj[j],  dtype= np.int)
        coarse_coef[i,local_id] = equiv_perm(1, 1)/25
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
P = spsolve(coarse_coef,coarse_q)
end = time.time()
print("This step lasted {0}s".format(end-start))

for cvolume,index in zip(M.coarse_volumes,range(len(P))):
    M.pressure[cvolume.volumes.global_id[cvolume.volumes.all]] = P[index]

print("Printing results")
M.core.print()
