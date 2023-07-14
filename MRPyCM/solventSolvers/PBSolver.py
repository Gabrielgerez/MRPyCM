import numpy as np
from vampyr import vampyr3d as vp
from . import GPESolver


class PBSolver(GPESolver):
    
    
    def __init__(self, rho, eps, kappa, Poisson_operator, Derivative_operator):
        super().__init__(rho, eps, Poisson_operator, Derivative_operator)
        self.k_sq = kappa
        
        
    def computeGamma(self, V_tot, epsilon):
        gamma = super().computeGamma(V_tot, epsilon)

        sinh = vp.FunctionMap(np.sinh, epsilon)
        k_sq_sinh_V_tot = self.k_sq * sinh(V_tot)
        k_sq_sinh_V_tot.crop(epsilon)
        return gamma - (1/(4*np.pi))*k_sq_sinh_V_tot