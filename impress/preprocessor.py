import numpy as np
import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh

start = time.time()

M = msh('mesh/20.h5m', dim = 3)

end = time.time()

print("The preprocessing step lasted {0}s".format(end-start))
