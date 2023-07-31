from . import ShiftFunction


class DHScreening(ShiftFunction):
    """Square of the Debye Huckel Screening parameter
    dependent on three-dimensional space
    """
    
    
    def __call__(self, r):
        C_eval = self.C(r)
        return (1.0-C_eval)*(self.outside**2)
    

