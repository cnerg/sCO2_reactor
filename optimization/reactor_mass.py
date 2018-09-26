# import Flow class
from ht_functions import Flow, oned_flow_modeling
# import physical properties
import physical_properties as pp

def reactor_mass(fuel, coolant, power, m_dot, T, P,
                 clad='Inconel-718', refl='Carbon',
                 cool_r=0.005, clad_t=0.00031, AR=1):
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
    flow_props = pp.FlowProperties(coolant, power, m_dot, T, P)
    # initialize reactor model
    rxtr = Flow(cool_r, clad_t, AR, power, fuel, coolant, clad, refl, flow_props)
    # perform 1D calculation
    oned_flow_modeling(rxtr)
    
    return rxtr.mass
