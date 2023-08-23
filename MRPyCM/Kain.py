import numpy as np
from vampyr import vampyr3d as vp

class Kain():
    
    instances = []
    
    def __init__(self, history):
        self.instance_index = len(Kain.instances)
        self.instances.append(self)
        self.history = history
        
        self.A = np.zeros((1, 1))
        self.b = np.zeros((1))
        self.c = np.zeros((1))
        
        self.f = []
        self.df = []
    
    
    @property
    def A(self):
        return self.instances[self.instance_index]._A
    

    @A.setter
    def A(self, A):
        self.instances[self.instance_index]._A = A
        
    
    @property
    def b(self):
        return self.instances[self.instance_index]._b
    

    @b.setter
    def b(self, b):
        self.instances[self.instance_index]._b = b
        
    
    @property
    def c(self):
        return self.instances[self.instance_index]._c
    

    @c.setter
    def c(self, c):
        self.instances[self.instance_index]._c = c
        
    
    @property
    def f(self):
        return self.instances[self.instance_index]._f
    
    @f.setter
    def f(self, f):
        self.instances[self.instance_index]._f = f
    
    
    @property
    def df(self):
        return self.instances[self.instance_index]._df
    
    @df.setter
    def df(self, df):
        self.instances[self.instance_index]._df = df
    
        
    def accelerate(self, func, dfunc , prec):
        
        self.f.append(func)
        self.df.append(dfunc)
        
        if (len(self.f) == 1) or (len(self.df) == 1):
            return func, dfunc
        
        if ((len(self.f) >= self.history) or (len(self.df) >= self.history)):
            self.f.pop(0)
            self.df.pop(0)

        self.setupLinearSystem()
        self.solveLinearSystem()
        kain_update = self.expandSolution(prec)
        
        return self.f[-1], kain_update
        

    def setupLinearSystem(self):
        nHistory = len(self.f) -1

        # Compute matrix A
        self.A = np.zeros((nHistory, nHistory))
        self.b = np.zeros(nHistory)
        phi_m = self.f[nHistory]
        fPhi_m = self.df[nHistory]
        
        for i in range(nHistory):
            phi_i = self.f[i]
            dPhi_im = phi_i - phi_m
            
            for j in range(nHistory):
                fPhi_j = self.df[j]
                dfPhi_jm = fPhi_j - fPhi_m
                self.A[i, j] -= vp.dot(dPhi_im, dfPhi_jm)
            self.b[i] += vp.dot(dPhi_im, fPhi_m)
        return   


    def solveLinearSystem(self):
        self.c = np.zeros(len(self.b))
        self.c = np.linalg.solve(self.A, self.b)
        return


    def expandSolution(self, prec):
        nHistory = len(self.f) -1

        phi_m = self.f[nHistory].deepCopy()
        fPhi_m = self.df[nHistory].deepCopy()
        
        for j in range(nHistory):
            fPhi_m += (self.f[j] + self.df[j] - phi_m - fPhi_m)*self.c[j]
        
        return fPhi_m.deepCopy()
        


        
        
            