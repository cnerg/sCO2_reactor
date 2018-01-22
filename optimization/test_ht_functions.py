import math
from physical_constants import *
import pytest
from ht_functions import FlowIteration

# parameters for test case
diameter = 0.01
radius = diameter / 2.0
PD = 2
c = 0.00031
L = 0.5
N = 6

def test_set_geom():
    """Test geometry initialization.
    """
    test = FlowIteration(diameter, PD, c, L)
    test.guess = N
    # expected values
    exp_De = diameter
    exp_A_flow = math.pi * radius**2
    exp_A_fuel = math.sqrt(3) * ((radius + c) * PD * 2)**2 / 2.0 -\
            (radius + c)**2 * math.pi
    # calculated values
    test.set_geom()
    # compare
    assert exp_De == test.D_e
    assert exp_A_flow == test.A_flow
    assert exp_A_fuel == test.A_fuel

def test_characterize_flow():
    """Test flow characterization.
    """
    test = FlowIteration(diameter, PD, c, L)
    test.guess = N
    # get geom
    test.set_geom()
    #expected values
    exp_G = m_dot / (test.A_flow * N)
    exp_v = exp_G / rho_cool
    exp_Re = rho_cool * exp_v * test.D_e / mu
    exp_Pr = Cp_cool * mu / k_cool
    exp_Nu = 0.023*math.pow(exp_Re,0.8)*math.pow(exp_Pr,0.4)
    exp_f = 0.184 / math.pow(exp_Re,0.2)
    exp_h = exp_Nu * k_cool / test.D_e
    # calculated values
    test.characterize_flow()
    # compare
    assert exp_G == test.G_dot
    assert exp_v == test.v
    assert exp_Re == test.Re
    assert exp_Nu == test.Nu
    assert exp_f == test.f
    assert exp_h == test.h_bar

def test_q():
    """Test q_per_channel calculation.
    """
    test = FlowIteration(diameter, PD, c, L)
    test.guess = N
    # get geom and flow conditions
    test.set_geom()
    test.characterize_flow()
    #expected value
    exp_q_per_channel = 22608.9
    test.get_q_per_channel()
    # compare
    assert abs(exp_q_per_channel - test.q_per_channel) < 1.0
    
def test_dp():
    """Test subchannel dp calculation.
    """
    test = FlowIteration(diameter, PD, c, L)
    test.guess = N
    test.set_geom()
    test.characterize_flow()
    test.calc_dp()
    # expected dp
    exp_dp = test.f * L * rho_cool * test.v ** 2 / (2*test.D_e)
    # compare
    assert exp_dp == test.dp
