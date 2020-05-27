from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling
import numpy as np

lu = ParallelLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
#lu = ParallelLocalUpscaling(mesh_file = None, dataset ='imex_datasets/spe_10_case_2/spe_10_case_2.dat')

kefx, kefy, kefz = [], [], []
kx, ky, kz = [], [], []
coarse_problems = []

for a, b in zip(lu.effective_permeability, lu.distribution):
    for c, d in zip(a, b):
        coarse_problems.append(int(d))
        kefx.append(c[0])
        kefy.append(c[1])
        kefz.append(c[2])

for problem in coarse_problems:
    index = coarse_problems.index(problem)
    x, y, z = kefx[index], kefy[index], kefz[index]
    volumes = lu.coarse.elements[problem].volumes.global_id[:]
    perm = lu.mesh.permeability[volumes]

    # X direction error
    permx = perm[:,0]
    errorx = ((x - permx)/x)*100

    # Y direction error
    permy = perm[:,1]
    errory = ((y - permy)/y)*100

    # Z direction error
    permz = perm[:,2]
    errorz = ((z - permz)/z)*100


print('\nMaximum error in x direction is {}%'.format(np.array(errorx).max()))
print('Maximum error in y direction is {}%'.format(np.array(errory).max()))
print('Maximum error in z direction is {}%'.format(np.array(errorz).max()))
