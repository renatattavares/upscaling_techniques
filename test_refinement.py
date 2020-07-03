import numpy as np
from imex_integration.mesh_constructor import MeshConstructor as mc

mesh1 = mc(np.array([1,1,1]), np.array([1,1,1]), 'mesh1.h5m')
mesh2 = mc(np.array([2,2,2]), np.array([0.5,0.5,0.5]), 'mesh2.h5m')
