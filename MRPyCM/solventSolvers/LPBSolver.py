import numpy as np

from . import PBSolver


class LPBSolver(PBSolver):
    def computeGamma(self, V_tot, epsilon):
        gamma = super(PBSolver, self).computeGamma(V_tot, epsilon)
        k_sq_V_tot = self.k_sq * V_tot
        k_sq_V_tot.crop(epsilon)
        return gamma - k_sq_V_tot
