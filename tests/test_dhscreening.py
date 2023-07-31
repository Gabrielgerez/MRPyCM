from math import isclose
import MRPyCM
from vampyr import vampyr3d as vp

C = MRPyCM.Cavity([[0.0, 0.0, 0.0]], [1.0], 0.2)
MRA = vp.MultiResolutionAnalysis(box=[-5, 5], order=5)
I = 0.1
kappa_out = MRPyCM.computeKappaOut(2.0, I)
k_sq = MRPyCM.DHScreening(C, inside=0.0, outside=kappa_out)

epsilon = 1.0e-3
P_eps = vp.ScalingProjector(mra=MRA, prec=epsilon)
k_sq_tree = P_eps(k_sq)


def test_kappa_out_val():
    assert isclose(kappa_out, 0.3446303846489457, rel_tol=1e-09, abs_tol=1.0e-09)


def test_inside_val():
    assert isclose(k_sq([0.0, 0.0, 0.0]), 0.0, rel_tol=1e-09, abs_tol=1.0e-09)
    
    
def test_outside_val():
    assert isclose(k_sq([2.0, 2.0, 2.0]), kappa_out**2, rel_tol=1e-09, abs_tol=1.0e-09)
    

def test_inside_val_tree():
    assert isclose(k_sq_tree([0.0, 0.0, 0.0]), 0.0, rel_tol=epsilon, abs_tol=epsilon)
    
    
def test_outside_val_tree():
    assert isclose(k_sq_tree([2.0, 2.0, 2.0]), kappa_out**2, rel_tol=epsilon, abs_tol=epsilon)
