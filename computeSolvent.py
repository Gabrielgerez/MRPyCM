import MRPyCM.soluteFunctions as sf
import MRPyCM.solventSolvers as ss
import MRPyCM.utilities as ut
import numpy as np
from vampyr import vampyr3d as vp

if __name__ == '__main__':
    import sys


# Global parameters
def main(*args, **kwargs):
    """
    input parameter dictionary:
    "order" : int
    "box" : list of floats
    "prec" : float
    "cav_coords" : list of lists of floats
    "radii" : list of floats
    "sigma" : float
    "eps_out" : float
    "perm_type" : string
    "solvent_type" : string
    "I" : float
    """
    keys =kwargs.keys()
    k = kwargs["order"] if ("order" in keys) else 5                            # Polynomial order
    L = kwargs["box"] if ("box" in keys) else  [-10,10]                        # Simulation box size
    epsilon = kwargs["prec"] if ("prec" in keys) else 1.0e-4                   # Relative precision
    
    # Define MRA and multiwavelet projector
    MRA = vp.MultiResolutionAnalysis(order=k, box=L)
    P_eps = vp.ScalingProjector(mra=MRA, prec=epsilon)
    D_abgv = vp.ABGVDerivative(mra=MRA, a=0.0, b=0.0)
    Poissop = vp.PoissonOperator(mra=MRA, prec=epsilon)

    # nuclear density and total molecular density to compute the vacuum potential
    #TODO: implement a molecular density compute function in utilities
    beta = 10000
    alpha = (beta / np.pi)**(3.0/2.0)
    charge = 1.0
    dens = P_eps(vp.GaussFunc(beta=beta, alpha=alpha*charge, position=[0.0, 0.0, 0.0], poly_exponent=[0,0,0]))

    # Solvent part

    cav_coords = kwargs["cav_coords"] if ("cav_coords" in keys) else [[0.0000000000,    0.0000000000,    0.000000000]] #centered in 
    radii = kwargs["radii"] if ("radii" in keys) else [3.7794522492515403]  
    width = kwargs["sigma"] if ("sigma" in keys) else 0.2 

    C = sf.Cavity(cav_coords, radii, width)
    
      
    eps_out = kwargs["eps_out"] if ("eps_out" in keys) else 2.0
    
    if ("linear" == kwargs["perm_type"].lower()):
        perm = sf.Linear(C, inside=1.0, outside=eps_out)
    else:
        perm = sf.Exponential(C, inside=1.0, outside=eps_out)
        
    perm_tree = P_eps(perm)
    
        
    if ("pb" == kwargs["solvent_type"].lower()):
        I = kwargs["I"] if ("I" in keys) else 0.1
        kappa_sq = sf.DHScreening(C, inside=0.0, outside=ut.computeKappaOut(eps_out, I))
        k_sq_tree = P_eps(kappa_sq)
        
        PB_SCRF = ss.PBESolver(dens, perm_tree, k_sq_tree, Poissop, D_abgv)
        reaction_op = ss.Reaction_operator(PB_SCRF)
        reaction_op.setup(epsilon)
        
        print("E_R: ", reaction_op.trace())
        
    elif ("lpb" == kwargs["solvent_type"].lower()):
        I = kwargs["I"] if ("I" in keys) else 0.1
        kappa_sq = sf.DHScreening(C, inside=0.0, outside=ut.compute_kappa_out(eps_out, I))
        k_sq_tree = P_eps(kappa_sq)
        
        LPBE_SCRF = ss.LPBESolver(dens, perm_tree,k_sq_tree, Poissop, D_abgv)
        reaction_op = ss.Reaction_operator(LPBE_SCRF)
        reaction_op.setup(epsilon)
        
        print("E_R: ", reaction_op.trace())
        
    else:
        SCRF = ss.GPESolver(dens, perm_tree, Poissop, D_abgv)
        reaction_op = ss.Reaction_operator(SCRF)
        reaction_op.setup(epsilon)
        
        print("E_R: ", reaction_op.trace())

    pass  


if __name__ == '__main__':
    arg_dict = eval(sys.argv[1])
    print(arg_dict)
    main(**arg_dict)