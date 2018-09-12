# import Flow class
from ht_functions import Flow, oned_flow_modeling
# import physical properties
import physical_properties as pp


def reactor_mass(fuel, coolant, power, m_dot, T, P, 
                 cool_r=0.005, clad_t=0.00031, AR=1):
    """Produce the mass of a valid reactor design given flow conditions and
    required thermal power output.

        Arguments:
        ----------
    """

    # get coolant properties
    flow_props = pp.FlowProperties(coolant, power, m_dot, T, P)
    print(flow_props.__dict__)
    # initialize reactor model
    rxtr = Flow(cool_r, clad_t, AR, power, fuel, coolant, flow_props)
    # perform calculation
    oned_flow_modeling(rxtr)
    
    return rxtr.mass

mass = reactor_mass('UW', 'H2O', 200000, 0.2, (1000, 1000), (17.5e6, 17.3e6))
print(mass)
