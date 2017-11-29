""" This module determines maxinum heat generation in a sCO2 reactor fuel pin as
a function of pitch to diameter ratio, and mass flow.

The following functions are contained in this module:

"""
import matplotlib.pyplot as plt
import math
import copy
import numpy as np

# define global variables
D = 0.01 # flow channel diameter [m]
PD = 1.05 # pitch to diameter ratio [-] 
L = 0.5 # fuel length [m]
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
    
    q_trip = q_bar / (2 * r_f ** 2 * L)

    t_clad = T_centerline - q_trip * r_f ** 2 / (4 * k_f)
     
    return q_bar, t_clad, q_trip


def get_mass_flux_triangle(PD, N_elements):
    """Calculate coolant mass flux for a triangle lattice.
    Arguments: PD (float) pitch-to-diameter ratio
               N_elements (int) number of fuel elements
               D (float): diameter of fuel pin
               m_dot (float): core mass flow rate.
    Returns: G_dot (float): coolant mass flux
    """
    r_f_c = r_f + c
    diam = r_f_c * 2
    pitch = PD * diam
    A_channel = (pitch ** 2 * math.sqrt(3) - math.pi * 2 * r_f_c ** 2) / 4
    flow_area = A_channel * N_elements * 2
    flow_perim = math.pi * diam * N_elements
    G_dot = m_dot / flow_area
    D_e = 4 * flow_area / flow_perim
    v = G_dot / rho_cool
    
    return G_dot, D_e, v


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
    q_bar, t_clad, q_trip = get_q_bar(r_f, 1.0, dt, 0)
    N_pins = math.ceil(Q_therm / q_bar)
    G_dot, D_e, v = get_mass_flux_triangle(D/2, N_pins)
    h_bar, Re = get_h_bar(G_dot, D_e)

    q_bar_new = get_q_bar(r_f, h_bar, dt, 1)[0]
    N_pins_new = math.ceil(Q_therm / q_bar_new)
    dP = calc_dp(D_e, v, Re)
    while N_pins != N_pins_new:
        N_pins = copy.deepcopy(N_pins_new)
        q_bar, t_clad, q_trip = get_q_bar(r_f, h_bar, dt, 1)
        N_pins_new = math.ceil(Q_therm / q_bar)
        G_dot, D_e, v = get_mass_flux_triangle(PD, N_pins)
        h_bar, Re = get_h_bar(G_dot, D_e)
    
    return [N_pins, h_bar, Re, dP, q_bar, t_clad, q_trip]

if __name__ == '__main__':
    data = {'N_pins' : [],
            'h_bar' : [],
            'Re' : [],
            'dP' : [],
            'q_bar' : [],
            't_clad' : [],
            'q_trip' : []}

    PDs = np.linspace(1.1, 1.5, 25)
    for PD_ratio in PDs:
        results = run_iter(PD_ratio)
        for idx, type in enumerate(data.keys()):
            data[type].append(results[idx])
    plt.plot(PDs, data['q_bar'])
    plt.show()
