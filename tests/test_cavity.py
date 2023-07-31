from math import isclose
import MRPyCM
from vampyr import vampyr3d as vp


C = MRPyCM.Cavity([[0.0, 0.0, 0.0]], [1.0], 0.2)
MRA = vp.MultiResolutionAnalysis(box=[-5, 5], order=5)
epsilon = 1.0e-3
P_eps = vp.ScalingProjector(mra=MRA, prec=epsilon)
C_tree = P_eps(C)

C_2 = MRPyCM.Cavity([[-0.35, 0.0, 0.0], [0.35, 0.0, 0.0]], [1.0, 1.0], 0.2)
C_2_tree = P_eps(C_2)

def test_inside_val():
    assert isclose(C([0.0, 0.0, 0.0]), 1.0, rel_tol=1e-09, abs_tol=1e-09)
    

def test_outside_val():
    assert isclose(C([2.0, 2.0, 2.0]), 0.0, rel_tol=1e-09, abs_tol=1e-09)

def test_cavity_volume():
    C_vol = C_tree.integrate()
    assert isclose(C_vol, 4.440117592926352, rel_tol=1e-04, abs_tol=0.0)

def test_inside_val_tree():
    assert isclose(C_tree([0.0, 0.0, 0.0]), 1.0, rel_tol=1e-04, abs_tol=0.0)
    

def test_outside_val_tree():
    assert isclose(C_tree([2.0, 2.0, 2.0]), 0.0, rel_tol=1e-04, abs_tol=0.0)
    
    
def test_cavity_2_volume():
    C_2_vol = C_2_tree.integrate()
    assert isclose(C_2_vol, 6.774723601250772, rel_tol=1e-04, abs_tol=0.0)