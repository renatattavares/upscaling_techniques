import numpy as np
import pandas as pd

# phi = pd.read_excel('phi.xlsx')
# phi_array = phi.values
# phi_array = np.reshape(phi_array, newshape = (224400, 5))
# empty_array = np.zeros((1,5))
# phi = pd.DataFrame({'1': phi_array[:,0], '2': phi_array[:,1], '3': phi_array[:,2], '4': phi_array[:,3], '5': phi_array[:,4]})
# np.savetxt('phi_final.txt', phi.values, fmt = '%f')
#
perm = pd.read_excel('permk.xlsx')
perm_array = perm.values
missing_line = np.array([5.48424, 8.57781, 7.16766, 9.19119, 3.21411, 1.60526])
perm_array = np.append(missing_line, perm_array)
perm_array = np.reshape(perm_array, newshape = (224400, 5))
perm = pd.DataFrame({'1': perm_array[:,0], '2': perm_array[:,1], '3': perm_array[:,2], '4': perm_array[:,3], '5': perm_array[:,4]})
np.savetxt('permk_final.txt', perm.values, fmt = '%f')
