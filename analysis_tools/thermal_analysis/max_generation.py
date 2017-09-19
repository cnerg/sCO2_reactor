""" This module determines maximum heat generation in a sCO2 reactor fuel pin as
a function of pitch to diameter ratio, and mass flow.

The following functions are contained in this module:

"""

import math

# define global variables

L = 50 # fuel length [cm]
k_c = 21 # clad conductivity [W/m-k]
k_f = 21 # fuel conductivity [W/m-k]
c = 0.0031 # clad thickness [m] 

def get_q_s(R, h, dt):
    """Calculate heat transfer from fuel element.
    Arguments: R (float): fuel radius
               h (float): heat transfer coefficient
    Returns: q_s (float): total heat transfer from fuel element.
    """
    
    A_r = 2 * math.pi * R * L
    A_m = (2 * math.pi * c * L) / math.log((R + c) / R)
    A_r_c = 2 * math.pi * (R + c) * L
    R_prime = (R / (2 * k_f * A_r)) + (c / (k_c * A_m)) + (1 / (h * (R + c)))

    q_s = dt / R_prime

    return q_s

def get_mass_flux(PD, N_elements, D, m_dot):
    """Calculate coolant mass flux.
    Arguments: PD (float) pitch-to-diameter ratio
               N_elements (int) number of fuel elements
               D (float): diameter of fuel pin
               m_dot (float): core mass flow rate.
    Returns: G_dot (float): coolant mass flux
    """

    flow_area = (PD * D) ** 2 * (((D ** 2) * math.pi) / 4) * N_elements
    flow_perim = math.pi * D * N_elements

    G_dot = m_dot / flow_area
    D_e = 4 * flow_area / flow_perim 

    return G_dot, D_e
    
def get_h_bar():
    """Calculate average heat transfer coefficient.
    Arguments: 
    Returns: h_bar (float) average heat transfer coefficient.
    """
    
    Re = G_dot  * D_e / mu
    Pr = C_p * mu / k_cool

    # equation is combination of 9-30, 9-31a in El-Wakil's Nuclear Heat
    # Transport. Valid for 1.1 <= PD <= 1.3
    Nu = (0.042 * PD - 0.024) * Re ** 0.8 * Pr ** 0.333
    h_bar = Nu * k / D_e

    return h_bar

def get_q_trip():



if __name__ == 'main':

    get flow area
    get mass flux
    get flow velocity
    get reynolds
    get nusselt 
    get h_bar
    get q_s
    get q'''
