from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

############### RUN LOCAL PROBLEMS ###############
#lp = LocalProblems(mesh_file = 'mesh/20.h5m', boundary_condition_type = 1)

############### RUN LOCAL UPSCALING ###############
lu = LocalUpscaling(mesh_file = 'mesh/20.h5m', boundary_condition_type = 1)

############### RUN IMPRESS ###############
#mesh = impress(mesh_file = 'mesh/20.h5m', dim = 3)
