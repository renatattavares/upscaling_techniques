from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from upscaling_procedures.local_upscaling import LocalUpscaling

lu = LocalUpscaling(preprocessor = impress, coarse_config = coarse_config, mesh_file = 'mesh/20.h5m', boundary_condition_type = 1)
