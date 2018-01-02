"""Physical Constants Used for 1D simulation
*** All values are bulk flow values averaged axially across the core ***
"""
m_dot = 0.75 # coolant flow [kg/s]
k_clad = 173e-3 # clad conductivity: W [kW/m-k]
#k_fuel = 21e-3 # fuel conductivity: UN [kW/m-k]
k_fuel = 42e-3 # fuel conductivity: UN [kW/m-k]
k_cool = 0.07531e-3 # coolant conductivity [kW/m-k]
mu = 0.00004306 # coolant viscosity [kg/m-s]
Cp_cool = 1.274 # coolant specific heat [kJ/kg-k]
rho_cool = 87.13 # coolant density [kg/m^3]
Q_therm = 131.2 # core thermal power [kW]
T_in = 962.9 # core inlet temp [K] (from power cycle model)
T_out = 1100 # core outlet temp [K] (from power cycle model)
T_bulk = T_in + (T_out - T_in) / 2 # bulk coolant temp. [K]
T_centerline = 1847.5 # centerline fuel temperature [K]
P_in = 1.791e7 # inlet pressure [Pa]
P_out = 1.742e7 # outlet pressure [Pa]
dp_allowed = abs(P_out - P_in)
