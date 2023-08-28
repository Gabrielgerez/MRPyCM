from vampyr import vampyr3d as vp
import numpy as np
from .. import utilities as ut
#from ..Kain import Kain


class GPESolver():
    
    instances = []
    
    def __init__(self, rho, eps, Poisson_operator, Derivative_operator, prec, maxiter=100, hist=0):

        self.instance_index = len(GPESolver.instances)
        GPESolver.instances.append(self)

        self.iterations = 0
        self.P = Poisson_operator
        self.D = Derivative_operator
        self.Density = rho
        self.Permittivity = eps
        self.max_iter = maxiter
        self.hist = hist  # not used right now
        
        self.rho_eff = (4*np.pi)*((self.Density * (self.Permittivity)**(-1))  - (self.Density))
        self.V_vac = self.P((4*np.pi)*(self.Density))
        
        self.V_R = self.initReactionPotential(prec)
        
    
    def initReactionPotential(self, prec):
        """
        we want to do a run of the scrf solver only with the vacuum potential 
        to get an initial guess of the reaction potential. 
        Then compute the very first reaction potential
        """    
        
        gamma_0 = self.computeGamma(self.V_vac, prec)
        V_R_0 = self.P(self.rho_eff + gamma_0)
        return V_R_0.deepCopy()
        
    
        
    def computeGamma(self, V_tot, epsilon):
        gamma =  vp.dot(vp.gradient(self.D,V_tot), vp.gradient(self.D, self.Permittivity)) * ( self.Permittivity**(-1))
        gamma = gamma.crop(epsilon)
        return gamma.deepCopy()

    
    def setup(self, prec):
        #kain = Kain(self.hist)
        # start loop
        update = 1.0
        print(f"Iter.{' '*2}Norm{' '*12}Update{' '*10}Energy (a.u.){' '*3}Energy update (a.u.)")
        print(f"{'-'*75}")
        E_r_old = 0.0
        for i in range(self.max_iter):
            V_tot = self.V_vac + self.V_R
                        
            # compute the surface charge distribution gamma= 1/4pi * \\nabla log(eps)\\cdot\\nabla V_tot
            gamma = self.computeGamma(V_tot, prec)

            # solve the generalized poisson equation for V_R 
            V_R_np1 = self.P((self.rho_eff) + (gamma))
            dV_R = V_R_np1 - self.V_R
            
            #hopefully do KAIN here
            # not yet implemented properly
            # if (self.hist > 0):
            #      self.V_R, dV_R = kain.accelerate(self.V_R, dV_R, prec)
            #hopefully do KAIN here
            
            self.V_R =  self.V_R + dV_R
            #check convergence
            update = dV_R.norm()
            E_r = self.computeEnergy()
            dE_r = E_r - E_r_old
            E_r_old = E_r
            print(f"{i}{' '*6}{self.V_R.norm():14.7e}  {update:14.7e}  {E_r:14.7e}  {dE_r:14.7e}") 
           
            if (update < prec):
                self.iterations = i
                break
        else:
            print("WARNING: GPESolver did not converge")
        
        
    def computeEnergy(self):
        return (1/2)*vp.dot(self.V_R, self.Density)


    @property 
    def V_R(self):
        return self.instances[self.instance_index]._V_R


    @V_R.setter
    def V_R(self, V_r):
        self.instances[self.instance_index]._V_R = V_r.deepCopy()  ## probably have this be a deepcopy
    
    
    @property
    def iterations(self):
        return self.instances[self.instance_index]._iterations


    @iterations.setter
    def iterations(self, i):
            self.instances[self.instance_index]._iterations = i



    
    
    