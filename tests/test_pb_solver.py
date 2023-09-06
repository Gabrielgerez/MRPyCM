from vampyr import vampyr3d as vp
import MRPyCM
import numpy as np

k = 5
L = [-10, 10]
epsilon = 1.0e-4

MRA = vp.MultiResolutionAnalysis(order=k, box=L)
P_eps = vp.ScalingProjector(mra=MRA, prec=epsilon)
P_eps_perm = vp.ScalingProjector(mra=MRA, prec=epsilon/100)
D_abgv = vp.ABGVDerivative(mra=MRA, a=0.0, b=0.0)
Poissop = vp.PoissonOperator(mra=MRA, prec=epsilon)

charge_coords = [[0.0, 0.0, 0.0]]
charges = [1.0]
charge_width = 1000

max_iter = 1

C = MRPyCM.Cavity(charge_coords, [1.0], 0.2)

dens = P_eps(MRPyCM.constructChargeDensity(charge_coords, charges, width_parameter=charge_width))

perm = P_eps_perm(MRPyCM.Exponential(C, inside=1.0, outside=2.0))
k_sq = P_eps_perm(MRPyCM.DHScreening(C, inside=0.0, outside=MRPyCM.computeKappaOut(eps_out=2.0, Ionic_strength=0.001)))

Solver = MRPyCM.PBSolver(dens, perm, k_sq, Poissop, D_abgv, epsilon, max_iter=max_iter)
    
def test_V_R_norm():
    assert np.isclose(Solver.V_R.norm(), 9.103739250292987, rtol=epsilon, atol=epsilon)

def test_V_R_integrate():
    assert np.isclose(Solver.V_R.integrate(), -729.339043117761, rtol=epsilon, atol=epsilon)
    
def test_zeroth_energy():
    assert np.isclose(Solver.computeEnergy(), -0.3607092953686604, rtol=epsilon, atol=epsilon)
    
def test_first_iteration():
    Solver.setup(epsilon)
    assert np.isclose(Solver.computeEnergy(), -0.2429892725523144, rtol=epsilon, atol=epsilon)

