from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from upscaling_procedures import local_upscaling

lu = local_upscaling('20.h5m', mesh = msh)
