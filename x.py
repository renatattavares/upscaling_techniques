# VECTOR FOR SETTING BOUNDARY CONDITIONS IN DIRECTON X
import numpy as np

elements = np.array((), dtype=int)

for i in range(25):
    elements = np.append(elements, ((i+1)-1)*5)
    print(((i+1)-1)*5)

elements2 = np.array((), dtype=int)

for i in range(25):
    elements2 = np.append(elements2, (4 + ((i+1)-1)*5))
    print(4 + ((i+1)-1)*5)
