import sys
sys.path.append('../optimization/')

# import Flow class
from ht_functions import Flow, oned_flow_modeling
# import physical properties
import physical_properties as pp
from mcnp_inputs import HomogeneousInput


def gen_reactor(fuel, coolant, power, m_dot, T, P,
                clad='Inconel-718', refl='Carbon',
                cool_r=0.0025, clad_t=0.00031, AR=1):
    """Produce the mass of a valid (coolable and critical) reactor design 
    given flow conditions and required thermal power output.

    Arguments:
    ----------
        fuel (str): fuel type
        coolant (str): coolant type
        power (float): thermal power [W]
        m_dot (float): mass flow rate [kg/s]
        T ((float, float)): inlet, outlet temp [K]
        P ((float, float)): inlet, outlet pressure [Pa]
        clad (optional, str): cladding type, default is Inconel-718
        refl (optional, str): reflector type, default is Carbon
        cool_r (optional, float): coolant channel radius [m]
        clad_t (optional, float): clad thickness [m]
        AR (optional, float): core aspect ratio [-]
    
    Returns:
    --------
        rxtr.mass (class atribute, float): reactor mass [kg]
    """
    # get coolant properties
    flow_props = pp.FlowProperties(coolant, m_dot, (T[0],T[1]), (P[0], P[1]))
    # initialize reactor model
    rxtr = Flow(cool_r, clad_t, AR, power, fuel, coolant, clad, refl, flow_props)
    # perform 1D calculation
    oned_flow_modeling(rxtr)

    return rxtr

def write_inp(rxtr):
    """
    """
    cool = rxtr.coolant
    fuel = rxtr.fuel
    ref_mult = rxtr.refl_m[fuel][cool] - 1
    matr = None
    if rxtr.fuel == 'UW':
        fuel = 'UN'
        matr = 'W'
    
    config = {'fuel' : fuel,
              'matr' : matr,
              'cool' : cool,
              'r_cool' : rxtr.r_channel*1e2,
              'clad' : 'Inconel-718',
              'rho_cool' : rxtr.fps.rho*1e-3,
              'fuel_frac' : rxtr.fuel_frac,
              'ref_mult' : ref_mult,
              'core_r' : 13.13, #rxtr.core_r*1e2,
              'power' : rxtr.gen_Q*1e-3
             }

    filename = '{0}-{1}.i'.format(fuel, cool)
    input = HomogeneousInput(config=config) 
    input.homog_core()
    input.write_input(filename)

SSnear = gen_reactor('UO2', 'CO2', 1.678e5, 1.2722, (793.8,900), (1.7906e7,1.7423e7))
write_inp(SSnear)
print(SSnear.__dict__)

SSnear = gen_reactor('UW', 'CO2', 1.678e5, 1.2722, (793.8,900), (1.7906e7,1.7423e7))
write_inp(SSnear)
print(SSnear.__dict__['core_r'])
print(SSnear.__dict__['fuel_frac'])
