""" This module determines maximum heat generation in a sCO2 reactor fuel pin as
a function of pitch to diameter ratio, and mass flow.

The following functions are contained in this module:

"""
import matplotlib.pyplot as plt
import math

# define global variables
#D = 0.1 # fuel diameter [cm]
PD = 1.1 # pitch to diameter ratio [-] 
L = 50 # fuel length [cm]
k_c = 21 # clad conductivity [W/m-k]
k_f = 21 # fuel conductivity [W/m-k]
k_cool = 0.07104 # coolant conductivity [W/m-k], choose min of inlet/outlet k
mu = 0.0000412 # coolant viscosity [kg/m-s], chose min of inlet/outlet
Cp_cool = 1.262 # coolant specific heat [kJ/kg-k], choose min of in/out
rho_cool = 80.52 # coolant density [kg/m^3], choose min of in/out
c = 0.0031 # clad thickness [m] 
Q_therm = 131200 # core thermal power [kW]
T_coolant = 1100 # bulk coolant temperature [K] (from power cylce model)
T_centerline = 1473.15 # centerline fuel temperature [K]
# ref: Mononitride Fuel for Fast Reactors (see Zotero)
m_dot = 0.75 # coolant mass flow rate [kg/s]

def get_q_s(R, h, dt):
    """Calculate heat transfer from fuel element.
    Arguments: R (float): fuel radius
               h (float): heat transfer coefficient
    Returns: q_s (float): total heat transfer from fuel element.
    """
    A_r = 2 * math.pi * R * L
    A_m = (2 * math.pi * c * L) / math.log((R + c) / R)
    A_r_c = 2 * math.pi * (R + c) * L
    R_prime = (R / (2 * k_f * A_r)) + (c / (k_c * A_m)) + (1 / (h * A_r_c))
    
    q_s = dt / R_prime
    
    return q_s

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
    flow_area = pitch_area - pin_area
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
    
    Re = rho_cool * D_e * v / mu
    Pr = Cp_cool * mu / k_cool
    
    # equation is combination of 9-30, 9-31a in El-Wakil's Nuclear Heat
    # Transport. Valid for 1.1 <= PD <= 1.3
    Nu = (0.042 * PD - 0.024) * math.pow(Re, 0.8) * math.pow(Pr, 0.333)
    h_bar = Nu * k_cool / D_e
    
    return h_bar

def get_q_trip(R, q_s):
    """Calculate volumetric heat generation.
    """
    
    q_trip = q_s / (R ** 2 * L * math.pi)

    return q_trip

def run_1D_analysis(n_pins):
    global D
    dt = T_centerline - T_coolant
    G_dot, D_e, v = get_mass_flux(PD, n_pins, D)
    h_bar = get_h_bar(G_dot, D_e, v)
    q_s = get_q_s(D / 2.0, h_bar, dt)
    q_trip = get_q_trip(D / 2.0, q_s)

    core_power = (D ** 2 * math.pi * L * q_trip / 4) * n_pins
    return core_power, q_trip, h_bar, q_s
    
def iterate_pins(n_pins_guess=1):
    
    core_power, q_trip, h_bar, q_s = run_1D_analysis(n_pins_guess)
    h_save = [h_bar]
    while core_power < Q_therm:
        n_pins_guess += 1
        core_power, q_trip, h_bar, q_s = run_1D_analysis(n_pins_guess)
        h_save.append(h_bar)
    
    #plt.plot(h_save)
    #plt.show()
    return n_pins_guess, q_trip, q_s

    

if __name__ == '__main__':

    global D
    q_s_save = []
    x_save = []
    for x in range(5, 100, 1):
        D = x / 100.0
        n_pins, q_trip, q_s = iterate_pins()
        q_s_save.append(q_s)
        x_save.append(D)
        print n_pins
    plt.plot(x_save, q_s_save)
    plt.show()
