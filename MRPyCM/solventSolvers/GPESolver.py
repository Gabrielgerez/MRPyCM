from vampyr import vampyr3d as vp
import numpy as np
from .. import utilities as ut


class GPESolver():
    
    instances = []
    
    def __init__(self, rho, eps, Poisson_operator, Derivative_operator, prec, maxiter=100):

        self.instance_index = len(GPESolver.instances)
        GPESolver.instances.append(self)

        self.iterations = 0
        self.P = Poisson_operator
        self.D = Derivative_operator
        self.Density = rho
        self.Permittivity = eps
        self.max_iter = maxiter
        
        self.rho_eff = -4*np.pi*(self.Density * (self.Permittivity)**(-1)  - self.Density)
        rho_vac = self.Density*(-4*np.pi)
        self.V_vac = self.P(rho_vac)
        
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
        gamma =  ((-1)*vp.dot( vp.gradient(self.D, self.Permittivity), vp.gradient(self.D,V_tot)) * ( self.Permittivity**(-1)))
        gamma = gamma.crop(epsilon)
        return gamma.deepCopy()
      
    
    @property 
    def V_R(self):
        return self._V_R


    @V_R.setter
    def V_R(self, V_r):
        self._V_R = V_r.deepCopy()  ## probably have this be a deepcopy

    
    def setup(self, prec):
        # start loop
        update = 1.0
        print("iter\t|V_R norm\t\t|update\t\t|E_r\t\t|E_r update") #TODO: make this a better print statement
        E_r_old = 0.0
        for i in range(self.max_iter):
            ut.plotFunction(self.V_R, left=-5.0, right=5.0, npoints=1000, title=f"V_R_{i}")
            
            V_tot = self.V_vac + self.V_R
            ut.plotFunction(V_tot, left=-5.0, right=5.0, npoints=1000, title=f"V_tot_{i}")
            
            # compute the surface charge distribution gamma= 1/4pi * \\nabla log(eps)\\cdot\\nabla V_tot
            gamma = self.computeGamma(V_tot, prec)

            # solve the generalized poisson equation for V_R 
            V_R_np1 = self.P((self.rho_eff) + (gamma))
            dV_R = V_R_np1 - self.V_R
            
            #hopefully do KAIN here
            #hopefully do KAIN here
            
            self.instances[self.instance_index].V_R += dV_R
            
            #check convergence
            update = dV_R.norm()
            E_r = self.computeEnergy()
            dE_r = E_r - E_r_old
            E_r_old = E_r
            print(f"{i} \t |{self.V_R.norm()}\t |{update} \t |{E_r} \t |{dE_r}") #TODO: make this a better print statement
           
            if (update < prec):
                self.instances[self.instance_index].iterations = i
                break
        else:
            print("WARNING: GPESolver did not converge")
        
        
    def computeEnergy(self):
        return (1/2)*vp.dot(self.V_R, self.Density)
        