from math import isclose
import MRPyCM
from vampyr import vampyr3d as vp

C = MRPyCM.Cavity([[0.0, 0.0, 0.0]], [1.0], 0.2)
MRA = vp.MultiResolutionAnalysis(box=[-5, 5], order=5)
perm_out = 2.0
lin_perm = MRPyCM.Linear(C, outside=perm_out)
exp_perm = MRPyCM.Exponential(C, outside=perm_out)
epsilon = 1.0e-3
P_eps = vp.ScalingProjector(mra=MRA, prec=epsilon)
lin_tree = P_eps(lin_perm)
exp_tree = P_eps(exp_perm)

def test_inside_val_linear():
    assert isclose(lin_perm([0.0, 0.0, 0.0]), 1.0, rel_tol=1e-09, abs_tol=0.0)
    

def test_inside_val_exponential():
    assert isclose(exp_perm([0.0, 0.0, 0.0]), 1.0, rel_tol=1e-09, abs_tol=0.0)
    
    
def test_outside_val_linear():
    assert isclose(lin_perm([2.0, 2.0, 2.0]), perm_out, rel_tol=1e-09, abs_tol=0.0)


def test_outside_val_exponential():
    assert isclose(exp_perm([2.0, 2.0, 2.0]), perm_out, rel_tol=1e-09, abs_tol=0.0)
    

def test_inside_val_linear_tree():
    assert isclose(lin_tree([0.0, 0.0, 0.0]), 1.0, rel_tol=epsilon, abs_tol=0.0)
    

def test_inside_val_exponential_tree():
    assert isclose(exp_tree([0.0, 0.0, 0.0]), 1.0, rel_tol=epsilon, abs_tol=0.0)
    
    
def test_outside_val_linear_tree():
    assert isclose(lin_tree([2.0, 2.0, 2.0]), perm_out, rel_tol=epsilon, abs_tol=0.0)


def test_outside_val_exponential_tree():
    assert isclose(exp_tree([2.0, 2.0, 2.0]), perm_out, rel_tol=epsilon, abs_tol=0.0)