""" This module determines maxinum heat generation in a sCO2 reactor fuel pin as
a function of pitch to diameter ratio, and mass flow.

The following functions are contained in this module:

"""
import matplotlib.pyplot as plt
import math

# define global variables
D = 0.01 # fuel diameter [m]
r_f = D / 2.0 # fuel radius [m]
PD = 1.1 # pitch to diameter ratio [-] 
L = 0.5 # fuel length [m]
k_c = 21e-3 # clad conductivity [kW/m-k]
k_f = 21e-3 # fuel conductivity [kW/m-k]
k_cool = 0.07527e-3 # coolant conductivity [kW/m-k], choose min of inlet/outlet k
nu = 0.00004305 # coolant viscosity [kg/m-s], chose min of inlet/outlet
Cp_cool = 1.274 # coolant specific heat [kJ/kg-k], choose min of in/out
rho_cool = 85.99 # coolant density [kg/m^3], choose min of in/out
c = 0.00031 # clad thickness [m] 
Q_therm = 131.2 # core thermal power [kW]
T_in = 962.9 # core inlet temp [K] (from power cycle model)
T_out = 1100 # core outlet temp [K] (from power cycle model)
T_bulk = T_in + (T_out - T_in) / 2 # bulk coolant temp. [K]
T_centerline = 1473.15 # centerline fuel temperature [K]
# ref: Mononitride Fuel for Fast Reactors (see Zotero)
m_dot = 0.75 # coolant mass flow rate [kg/s]

def get_q_bar(R, h, dt, h_bar_switch=1):
    """Calculate heat transfer from fuel element.
    Arguments: R (float): fuel radius
               h (float): heat transfer coefficient
    Returns: q_s (float): total heat transfer from fuel element.
    """
    
    R_prime = (1/k_f) + (1.0 / (2 * k_c)) * math.log((r_f + c) / r_f ) +\
            h_bar_switch / (2 * h * (r_f + c))

    q_bar = (dt * 2*L) / R_prime
    
    return q_bar

def get_mass_flux(PD, N_elements, D):
    """Calculate coolant mass flux.
    Arguments: PD (float) pitch-to-diameter ratio
               N_elements (int) number of fuel elements
               D (float): diameter of fuel pin
               m_dot (float): core mass flow rate.
    Returns: G_dot (float): coolant mass flux
    """
    pitch_area = (PD * D) ** 2.0 
    pin_area = math.pi * (D / 2.0) ** 2.0
    flow_area = ((pitch_area - pin_area) / 1e4) * N_elements
    flow_perim = math.pi * D * N_elements
    G_dot = m_dot / flow_area
    D_e = 4 * flow_area / flow_perim
    v = G_dot / rho_cool

    return G_dot, D_e, v
    
def get_h_bar(G_dot, D_e, v):
    """Calculate average heat transfer coefficient.
    Arguments: 
    Returns: h_bar (float) average heat transfer coefficient.
    """
    
    Re = rho_cool * D_e * v / nu
    Pr = Cp_cool * nu / k_cool
    
    # equation is combination of 9-30, 9-31a in El-Wakil's Nuclear Heat
    # Transport. Valid for 1.1 <= PD <= 1.3
    Nu = (0.042 * PD - 0.024) * math.pow(Re, 0.8) * math.pow(Pr, 0.333)
    h_bar = Nu * k_cool / D_e
    
    return h_bar

if __name__ == '__main__':
    dt = T_centerline - T_bulk
    # get first estimate for q_bar
    q_bar = get_q_bar(D / 2.0, 1.0, dt, 0)
    print q_bar
    N_pins = math.ceil(Q_therm / q_bar)
    print N_pins

