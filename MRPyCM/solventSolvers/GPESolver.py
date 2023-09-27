import sys
from vampyr import vampyr3d as vp
import numpy as np
from .. import utilities as ut

sys.path.insert(1, 'external/response/src/')
from KAIN import KAIN




class GPESolver():
    V_R : vp.FunctionTree
    iterations : int
    
    def __init__(self, rho, eps, Poisson_operator, Derivative_operator, maxiter=100, hist=0):

        self.Density = rho
        self.Permittivity = eps
        self.P = Poisson_operator
        self.D = Derivative_operator
        self.max_iter = maxiter
        self.hist = hist  # not used right now
        
        self.rho_eff = (4*np.pi)*((self.Density * (self.Permittivity)**(-1))  - (self.Density))
        self.iterations = 0
       
        self.V_vac = rho ## dummy value
        self.V_R = rho ## dummy value
            
    
    def initReactionPotential(self, prec):
        """
        we want to do a run of the scrf solver only with the vacuum potential 
        to get an initial guess of the reaction potential. 
        Then compute the very first reaction potential
        """    
        
        gamma_0 = self.computeGamma(self.V_vac, prec)
        self.V_R = self.P(self.rho_eff + gamma_0)
    
        
    def computeGamma(self, V_tot, epsilon):
        gamma =  vp.dot(vp.gradient(self.D,V_tot), vp.gradient(self.D, self.Permittivity)) * ( self.Permittivity**(-1))
        gamma = gamma.crop(epsilon)
        return gamma.deepCopy()

    
    def solveEquation(self, prec):
        print(f"Iter.{' '*2}Norm{' '*12}Update{' '*10}Energy (a.u.){' '*3}Energy update (a.u.)")
        print(f"{'-'*75}")
        
        self.V_vac = self.P((4*np.pi)*(self.Density))
        self.initReactionPotential(prec)
        
        Kain = KAIN(self.hist)
        # start loop
        update = self.V_R.norm()
        
        
        E_r = self.computeEnergy()
        dE_r = E_r
        print(f"{0:2d}{' '*5}{self.V_R.norm():14.7e}  {update:14.7e}  {E_r:14.7e}  {dE_r:14.7e}") 
        for i in range(self.max_iter):
            V_tot = self.V_vac + self.V_R
                        
            # compute the surface charge distribution gamma= 1/4pi * \\nabla log(eps)\\cdot\\nabla V_tot
            gamma = self.computeGamma(V_tot, prec)

            # solve the generalized poisson equation for V_R 
            V_R_np1 = self.P((self.rho_eff) + (gamma))
            dpotential = [V_R_np1 - self.V_R]
            #hopefully do KAIN here
            # not yet implemented properly
            if (self.hist > 0):
                dpotential = Kain.accelerate([self.V_R], dpotential)
            #hopefully do KAIN here
            
            self.V_R =  self.V_R + dpotential[0]
            #check convergence
            update = dpotential[0].norm()
            E_r_np1 = self.computeEnergy()
            
            dE_r = E_r_np1 - E_r
            E_r = E_r_np1
            print(f"{i+1:2d}{' '*5}{self.V_R.norm():14.7e}  {update:14.7e}  {E_r:14.7e}  {dE_r:14.7e}") 
            
            if (update < prec):
                self.iterations = i+1
                break
        else:
            print("WARNING: GPESolver did not converge")
        
        
    def computeEnergy(self):
        return (1/2)*vp.dot(self.V_R, self.Density)


    
    
    