from vampyr import vampyr3d as vp
import MRPyCM
import numpy as np


if __name__ == '__main__':
    import sys


# Global parameters
def run(*args, **kwargs):
    """
    input parameter dictionary:
    "order" : int
    "box" : list of floats
    "prec" : float
    "charge_width" : float
    "charges" : list of floats
    "charge_coords" : list of lists of floats
    "cav_coords" : list of lists of floats
    "cav_radii" : list of floats
    "boundary_width" : float
    "eps_out" : float
    "perm_formulation" : string
    "solvent_type" : string
    "ionic_strength" : float
    "kain_hist" : int
    "max_iter" : int
    """
    # Define parameters and defaults
    keys =kwargs.keys()
    k = kwargs["order"] if ("order" in keys) else 5                            # Polynomial order
    L = kwargs["box"] if ("box" in keys) else  [-10,10]                        # Simulation box size
    epsilon = kwargs["prec"] if ("prec" in keys) else 1.0e-4                   # Relative precision
    
    charge_width = kwargs["charge_width"] if ("charge_width" in keys) else 1000.0
    charges = kwargs["charges"] if ("charges" in keys) else [1.0]
    charge_coords = kwargs["charge_coords"] if ("charge_coords" in keys) else [[0.0000000000,    0.0000000000,    0.000000000]]
    
    cav_coords = kwargs["cav_coords"] if ("cav_coords" in keys) else  charge_coords 
    cav_radii = kwargs["cav_radii"] if ("cav_radii" in keys) else [1.0]  
    boundary_width = kwargs["boundary_width"] if ("boundary_width" in keys) else 0.2 
    
    eps_out = kwargs["eps_out"] if ("eps_out" in keys) else 2.0
    perm_formulation = kwargs["perm_formulation"] if "perm_formulation" in keys else "exponential"
    
    solvent_type = kwargs["solvent_type"] if ("solvent_type" in keys) else "gpe"
    ionic_strength = kwargs["ionic_strength"] if ("ionic_strength" in keys) else 0.1
    
    max_iter = kwargs["max_iter"] if ("max_iter" in keys) else 100
    kain_hist = kwargs["kain_hist"] if ("kain_hist" in keys) else 0
    
    
    # Define MRA and multiwavelet projector
    MRA = vp.MultiResolutionAnalysis(order=k, box=L)
    P_eps = vp.ScalingProjector(mra=MRA, prec=epsilon)
    D_abgv = vp.ABGVDerivative(mra=MRA, a=0.0, b=0.0)
    Poissop = vp.PoissonOperator(mra=MRA, prec=epsilon)

    # nuclear density and total molecular density to compute the vacuum potential
    dens = P_eps(MRPyCM.constructChargeDensity(charge_coords, charges, width_parameter=charge_width))

    # Solvent part
    C = MRPyCM.Cavity(cav_coords, cav_radii, boundary_width)
    
    if ("linear" == perm_formulation.lower()):
        perm = P_eps(MRPyCM.Linear(C, inside=1.0, outside=eps_out))
    else:
        perm = P_eps(MRPyCM.Exponential(C, inside=1.0, outside=eps_out))
        
        
    if ("pb" == solvent_type.lower()):
        k_sq = P_eps(MRPyCM.DHScreening(C, inside=0.0, outside=MRPyCM.computeKappaOut(eps_out, ionic_strength)))
        Solver = MRPyCM.PBSolver(dens, perm, k_sq, Poissop, D_abgv, epsilon, max_iter=max_iter, hist=kain_hist)
        
    elif ("lpb" == solvent_type.lower()):
        k_sq = P_eps(MRPyCM.DHScreening(C, inside=0.0, outside=MRPyCM.computeKappaOut(eps_out, ionic_strength)))
        Solver = MRPyCM.LPBSolver(dens, perm, k_sq, Poissop, D_abgv, epsilon, max_iter=max_iter, hist=kain_hist)
        
    else:
        Solver = MRPyCM.GPESolver(dens, perm, Poissop, D_abgv, epsilon, maxiter=max_iter, hist=kain_hist)
    
    
    reaction_op = MRPyCM.ReactionOperator(Solver)
    reaction_op.setup(epsilon)
    E_R = reaction_op.trace()
    print("E_R: ", E_R)
    return E_R, Solver.iterations
    


if __name__ == '__main__':
    arg_dict = eval(sys.argv[1])
    run(**arg_dict)