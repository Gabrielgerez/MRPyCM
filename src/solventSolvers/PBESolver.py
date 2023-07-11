import numpy as np
from vampyr import vampyr3d as vp
from . import GPESolver


class PBESolver(GPESolver):
    
    
    def __init__(self, rho, eps, kappa, Poisson_operator, Derivative_operator):
        super().__init__(rho, eps, Poisson_operator, Derivative_operator)
        self.k_sq = kappa
        
        
    def computeGamma(self, V_tot, epsilon):
        gamma = super().computeGamma(V_tot)

        sinh = vp.FunctionMap(np.sinh, epsilon)
        
        return gamma - (1/(4*np.pi))*self.k_sq * sinh(V_tot)