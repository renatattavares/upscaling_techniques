from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from upscaling_procedures.local_problems import LocalProblems

lp = LocalProblems(preprocessor = impress, coarse_config = coarse_config, mesh_file = 'mesh/25.h5m', boundary_condition_type = 1)
