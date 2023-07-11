import numpy as np 
from . import PBESolver


class LPBESolver(PBESolver):
    
    
    def computeGamma(self, V_tot, epsilon):
        gamma =  super(PBESolver, self).computeGamma(V_tot)
        return gamma - (1/(4*np.pi))*(self.k_sq * V_tot)