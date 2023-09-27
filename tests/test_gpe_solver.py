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
kain_hist = 0

C = MRPyCM.Cavity(charge_coords, [1.0], 0.2)

dens = P_eps(MRPyCM.constructChargeDensity(charge_coords, charges, width_parameter=charge_width))

perm = P_eps_perm(MRPyCM.Exponential(C, inside=1.0, outside=2.0))

Solver = MRPyCM.GPESolver(dens, perm, Poissop, D_abgv, max_iter=max_iter, hist=kain_hist)
Solver.iterateEquation(epsilon)

def test_V_vac_norm():
    assert np.isclose(Solver.V_vac.norm(), 12.375991947513185, rtol=epsilon, atol=epsilon)

def test_V_vac_integrate():
    assert np.isclose(Solver.V_vac.integrate(), 952.0261858394258, rtol=epsilon, atol=epsilon)
    
def test_V_R_norm():
    assert np.isclose(Solver.V_R.norm(), 5.461235106973656, rtol=epsilon, atol=epsilon)

def test_V_R_integrate():
    assert np.isclose(Solver.V_R.integrate(), -430.30804111722665, rtol=epsilon, atol=epsilon)
    
def test_zeroth_energy():
    assert np.isclose(Solver.energy_convergence[0], -0.35397429984367895, rtol=epsilon, atol=epsilon)
    
def test_first_iteration():
    assert np.isclose(Solver.energy_convergence[-1], -0.24140649819582582, rtol=epsilon, atol=epsilon)

