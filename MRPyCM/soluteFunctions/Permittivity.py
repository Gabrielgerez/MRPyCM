import numpy as np
from . import ShiftFunction

class Linear(ShiftFunction):
    
    
    def __call__(self, r):
        C_eval = self.C(r)
        permittivity = C_eval*(self.inside - self.outside) + self.outside
        return permittivity


class Exponential(ShiftFunction):
    
    
    def __call__(self, r):
        C_eval = self.C(r)
        permittivity = self.inside*np.exp((np.log((self.inside/self.outside)))*(1.0 - C_eval))
    
        return permittivity