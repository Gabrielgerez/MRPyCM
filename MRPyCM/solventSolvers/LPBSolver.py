import numpy as np 
from . import PBSolver


class LPBSolver(PBSolver):
    
    
    def computeGamma(self, V_tot, epsilon):
        gamma =  super(PBSolver, self).computeGamma(V_tot, epsilon)
        return gamma - (1/(4*np.pi))*(self.k_sq * V_tot)