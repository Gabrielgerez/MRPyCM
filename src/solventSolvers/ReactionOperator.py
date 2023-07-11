class ReactionOperator():
    
    
    def __init__(self, Solver):
        print("Initializing reaction operator")
        self.Solver = Solver
    
        
    def setup(self, prec):
        print("Setting up reaction operator")
        self.Solver.setup(prec)


    def trace(self):
        return self.Solver.computeEnergy()