from vampyr import vampyr3d as vp
import numpy as np


class GPESolver():
    
    instances = []
    
    
    def __init__(self, rho, eps, Poisson_operator, Derivative_operator):
        self.Density = rho
        self.Permittivity = eps
        self.P = Poisson_operator
        self.D = Derivative_operator
        self.V_R = rho
        self.instance_index = len(GPESolver.instances)
        GPESolver.instances.append(self)
        
        
    def computeGamma(self, V_tot, epsilon):
        return (1.0/(4*np.pi)) * vp.dot( vp.gradient(self.D, self.Permittivity), vp.gradient(self.D,V_tot)) * ( self.Permittivity**(-1))
    
    
    @property
    def Density(self):
        return self._Density 
    
    
    @Density.setter
    def Density(self, rho):
        self._Density = rho
    
    
    @property 
    def V_R(self):
        return self._V_R


    @V_R.setter
    def V_R(self, V_r):
        self._V_R = V_r

    
    def setup(self, prec):
        # compute rho_eff = rho/eps

        rho_eff = self.Density * (self.Permittivity)**(-1)
        
        # compute vacuum potential V_vac = \\int rho(r')/|r-r'| dr'
        V_vac = self.P(self.Density)
        
        # start loop
        update = 1.0
        print("iter\t|V_R norm\t\t|update\t\t|E_r\t\t|E_r update") #TODO: make this a better print statement
        E_r_old = 0.0
        for i in range(100):
            
            if (i==0):
                V_tot = V_vac
            else:
                V_tot = V_vac + self.instances[self.instance_index].V_R
            
            # compute the surface charge distribution gamma= 1/4pi * \\nabla log(eps)\\cdot\\nabla V_tot
            # in the first iteration, V_tot = V_vac else V_tot = V_vac + V_R
            gamma = self.computeGamma(V_tot, prec)

            # solve the generalized poisson equation for V_tot
            V_tot_np1 = self.P(rho_eff + gamma)
            #substract V_vac from V_tot to get V_R
            V_R_np1 = V_tot_np1 - V_vac
            if (i==0):
                dV_R = V_R_np1
                self.instances[self.instance_index].V_R = V_R_np1
            else:
                dV_R = V_R_np1 - self.instances[self.instance_index].V_R
                self.instances[self.instance_index].V_R += dV_R
            #hopefully do KAIN here
            
            #check convergence
            update = dV_R.norm()
            E_r = self.computeEnergy()
            dE_r = E_r - E_r_old
            E_r_old = E_r
            print(f"{i} \t |{self.instances[self.instance_index].V_R.norm()}\t |{update} \t |{E_r} \t |{dE_r}") #TODO: make this a better print statement
           
            if (update < prec):
                break
        else:
            raise RuntimeError("Did not converge")
        
        
    def computeEnergy(self):
        return (1/2)*vp.dot(self.Density, self.instances[self.instance_index].V_R)
        