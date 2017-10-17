""" This module determines maxinum heat generation in a sCO2 reactor fuel pin as
a function of pitch to diameter ratio, and mass flow.

The following functions are contained in this module:

"""
import matplotlib.pyplot as plt
import math
import copy
import numpy as np

# define global variables
D = 0.01 # fuel diameter [m]
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

def get_mass_flux_square(PD, N_elements):
    """Calculate coolant mass flux for a square lattice.
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

def get_mass_flux_triangle(PD, N_elements):
    """Calculate coolant mass flux for a triangle lattice.
    Arguments: PD (float) pitch-to-diameter ratio
               N_elements (int) number of fuel elements
               D (float): diameter of fuel pin
               m_dot (float): core mass flow rate.
    Returns: G_dot (float): coolant mass flux
    """
    pitch_area = (1 / 4.0) * (PD * D) ** 2 * math.sqrt(3)
    pin_area = math.pi * (r_f) ** 2.0
    flow_area = (pitch_area - 0.5 * pin_area) * N_elements
    flow_perim = math.pi * D * N_elements
    G_dot = m_dot / flow_area
    D_e = 4 * flow_area / flow_perim
    v = G_dot / rho_cool
    return G_dot, D_e, v

def get_h_bar(G_dot, D_e):
    """Calculate average heat transfer coefficient.
    Arguments: 
    Returns: h_bar (float) average heat transfer coefficient.
    """
    
    Re = D_e * G_dot/ nu
    Pr = Cp_cool * nu / k_cool
    
    # equation is combination of 9-30, 9-31b in El-Wakil's Nuclear Heat
    # Transport. Valid for 1.1 <= PD <= 1.3
    Nu = (0.026 * PD - 0.006) * math.pow(Re, 0.8) * math.pow(Pr, 0.333)
    h_bar = Nu * k_cool / D_e
    
    return h_bar, Re

def calc_dp(D_e, v, Re):
    """Calculate pressure drop in triangle-pitch subchannel
    Arguments: D_e (float) [m] hydraulic diameter
               v   (float) [m/s] bulk flow velocity
    Returns:
               dp (float) [Pa] core pressure drop
    """
    # friction factor correlation (Todreas & Kazimi eq. 9109b)
    c = 0.1458 + 0.03632 * (PD - 1) - 0.03333 * (PD -1) ** 2
    f_s = c / math.pow(Re, 0.18)

    return f_s * L * rho_cool * v ** 2 / (D_e * 2)

def run_iter(PD, diam=D):
    global r_f
    r_f = diam / 2.0 # fuel radius [m]
    
    dt = T_centerline - T_bulk
    # get first estimate for q_bar
    q_bar = get_q_bar(r_f, 1.0, dt, 0)
    N_pins = math.ceil(Q_therm / q_bar)
    G_dot, D_e, v = get_mass_flux_triangle(PD, N_pins)
    h_bar, Re = get_h_bar(G_dot, D_e)

    q_bar_new = get_q_bar(r_f, h_bar, dt, 1)
    N_pins_new = math.ceil(Q_therm / q_bar_new)
    dP = calc_dp(D_e, v, Re)
    while N_pins != N_pins_new:
        N_pins = copy.deepcopy(N_pins_new)
        q_bar = get_q_bar(r_f, h_bar, dt, 1)
        N_pins_new = math.ceil(Q_therm / q_bar)
        G_dot, D_e, v = get_mass_flux_triangle(PD, N_pins)
        h_bar, Re = get_h_bar(G_dot, D_e)

    return [N_pins, h_bar, Re, dP, q_bar]

if __name__ == '__main__':
    data = {'N_pins' : [],
            'h_bar' : [],
            'Re' : [],
            'dP' : [],
            'q_bar' : []}

    PDs = np.linspace(1.1, 1.5, 25)
    for PD_ratio in PDs:
        results = run_iter(PD_ratio)
        for idx, type in enumerate(data.keys()):
            data[type].append(results[idx])
    plt.plot(PDs, data['dP'])
    plt.show()
