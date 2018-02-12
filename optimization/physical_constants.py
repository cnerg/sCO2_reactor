"""Physical Constants Used for 1D simulation
*** All values are bulk flow values averaged axially across the core ***
"""
import math


def fuel_cond(T):
    """Estimate CERMET fuel conductivity based on T. Use a correlation from Webb
    and Charit "Analytical Determination of thermal conductivity of W-UO2 and
    W-UN CERMET nuclear fuels. Correlation provided for 60% UN
    """

    kc = 1.841e-19*math.pow(T, 6) - 2.097e-15*math.pow(T, 5) +\
        9.721e-12*math.pow(T, 4) - 2.369e-8*math.pow(T, 3) +\
        3.283e-5*math.pow(T, 2) - 0.0267*T + 63.18

    return kc


###############################################################################
#                                                                             #
#                     Given Parameters From Power Cycle                       #
#                                                                             #
###############################################################################
m_dot = 0.75  # coolant flow [kg/s]
Q_therm = 131000  # core thermal power [W]
T_in = 962.9  # core inlet temp [K] (from power cycle model)
T_out = 1100  # core outlet temp [K] (from power cycle model)
T_bulk = T_in + (T_out - T_in) / 2  # bulk coolant temp. [K]
rho_cool = 87.13  # coolant density [kg/m^3]
Cp_cool = 1274  # coolant specific heat [J/kg-k]
mu = 0.00004306  # coolant viscosity [kg/m-s]
P_in = 1.79064e7  # inlet pressure [Pa]
P_out = 1.74229e7  # outlet pressure [Pa]

###############################################################################
#                                                                             #
#                            Literature Values                                #
#                                                                             #
###############################################################################
T_fuel_max = 1847.5  # centerline fuel temperature [K]
T_centerline = T_fuel_max
k_clad = 108.3  # clad conductivity: W @ 1473 K [W/m-K]
k_cool = 0.07531  # coolant conductivity [W/m-k]
# conservative estimate for thermal conductivity at fuel centerline temperature.
k_fuel = fuel_cond(T_centerline)
rho_W = 19250  # clad density [kg/m^3]
rho_UN = 11300  # fuel density [kg/m^3]
fuel_frac = 0.6  # fraction of fuel in CERMET

###############################################################################
#                                                                             #
#                   Calculated Parameters From Power Cycle                    #
#                                                                             #
###############################################################################
dp_allowed = P_in - P_out  # pressure drop limit [Pa]
Pr = Cp_cool * mu / k_cool  # coolant Prandtl numbe[-]
# mixed density for CERMET fuel
rho_fuel = fuel_frac*rho_UN + (1-fuel_frac)*rho_W
