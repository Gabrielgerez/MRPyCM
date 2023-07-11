import numpy as np
from scipy.special import erf

class Cavity(object):
    
    
    def __init__(self, cav_coords, radii, width):
        self.r_list = cav_coords # list of cavity centers. Normally, but not always, nucleus coordinates. 
        self.R_list = radii  # list of cavity radii. Normally, but not always, nucleus radii.
        self.sigma = width  # width of the cavity boundary
   

    def __call__(self,  r):
        """
        r is a list of floats of length 3, can be a numpy array as well
        """
        r_vec = np.array(r)
        C = 1.0
        for i, r_i in enumerate(self.r_list):
            r_vec_i = np.array(r_i)
            s_i = np.linalg.norm(r_vec_i - r_vec) - self.R_list[i]
            O_i = (1.0/2.0)*(1 + erf(s_i/self.sigma))
            C_i = 1 - O_i
            C *= 1 - C_i
        C = 1.0 - C
        return C