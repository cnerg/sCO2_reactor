import nose
import math
from physical_constants import *
from nose.tools import assert_equal, assert_false, assert_true
from ht_functions import FlowIteration

def test_check_converge():
    """Test the function to test for convergence.
    """
    test = FlowIteration(1, 2, 0.00031, 1, 1)
    assert_false(test.check_converge())
    test.N_channels = test.guess
    assert_true(test.check_converge())
    test.dp = 1e6
    assert_false(test.check_converge())

def test_mass_flux():
    """Test mass flux function.
    """
    test = FlowIteration(1, 2, 0.00031, 1, 1)
    exp_G = m_dot / (math.pi * 0.5**2)
    exp_De = 2*test.r_channel
    exp_v = exp_G / rho_cool

    test.mass_flux_channel()
    assert_equal(exp_De, test.D_e)
    assert_equal(exp_G, test.G_dot)
    assert_equal(exp_v, test.v)

def test_nondim():
    """Test non-dimensional flow calculations.
    """
    test = FlowIteration(1, 2, 0.00031, 1, 1)
    # get mass flux data
    test.mass_flux_channel()
    exp_Re = rho_cool * test.v * test.D_e / mu
    exp_Pr = Cp_cool * mu / k_cool
    exp_Nu = 0.023*math.pow(exp_Re,0.8)*math.pow(exp_Pr, 0.4)
    test.calc_nondim()
    assert_equal(exp_Re, test.Re)
    assert_equal(exp_Pr, test.Pr)
    assert_equal(exp_Nu, test.Nu)

def test_h_bar():
    """Test function to calculate heat transfer coefficient.
    """ 
    test = FlowIteration(1, 2, 0.00031, 1, 1)
    test.mass_flux_channel()
    test.calc_nondim()
    exp_hbar = test.Nu * k_cool / test.D_e

    test.get_h_bar()
    assert_equal(exp_hbar, test.h_bar)

"""
SOME FUNCTION TO TEST QBAR GOES HERE EVENTUALLY DO IT ALEX!!!!!!!
"""
