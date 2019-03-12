import numpy as np

class FlowBasedPermeability():
    def __init__(self, M, flux_direction, bc):
        area = 1
        total_flow = 0.0
        flow_rate = 0.0
        self.M = M
        self.bc = bc
        
        if flux_direction == 'x':
            effective_permeability = self.x()
        elif flux_direction == 'y':
            pass
        elif flux_direction == 'z':
            effective_permeability = self.z()

    def x(self):

        for v in range(len(bc.elements)):
             flow_rate =  + equiv_perm(perm[v], perm[v+1])*area*(M.coarse_volumes[i].pressure_coarse[np.array(v)]-M.coarse_volumes[i].pressure_coarse[np.array(v+1)])
             total_flow = total_flow + flow_rate

        effective_permeability = total_flow/((area*rx*ry)*(M.coarse_volumes[i].pressure_coarse[0]-M.coarse_volumes[i].pressure_coarse[1]))

        return effective_permeability

    def y(self):
        pass

    def z(self):

        for v in range(len(bc.elements)):
             flow_rate =  + equiv_perm(perm[v], perm[v+rx*ry])*area*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+rx*ry])
             total_flow = total_flow + flow_rate

        effective_permeability = total_flow/((area*rx*ry)*(M.coarse_volumes[i].pressure_coarse[v]-M.coarse_volumes[i].pressure_coarse[v+rx*ry]))

        return effective_permeability
