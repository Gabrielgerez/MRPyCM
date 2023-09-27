import numpy as np

from .shiftfunction import ShiftFunction


class LinPerm(ShiftFunction):
    def __call__(self, r):
        C_eval = self.C(r)
        permittivity = self.inside * C_eval + self.outside * (1.0 - C_eval)
        return permittivity


class ExpPerm(ShiftFunction):
    def __call__(self, r):
        C_eval = self.C(r)
        permittivity = self.inside * np.exp(
            (np.log((self.outside / self.inside))) * (1.0 - C_eval)
        )

        return permittivity
