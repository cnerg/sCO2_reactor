m_dot = 0.75
k_c = 21e-3 # clad conductivity [kW/m-k]
k_f = 21e-3 # fuel conductivity [kW/m-k]
k_cool = 0.07527e-3 # coolant conductivity [kW/m-k], choose min of inlet/outlet k
nu = 0.00004305 # coolant viscosity [kg/m-s], chose min of inlet/outlet
Cp_cool = 1.274 # coolant specific heat [kJ/kg-k], choose min of in/out
rho_cool = 85.99 # coolant density [kg/m^3], choose min of in/out


Q_therm = 131.2 # core thermal power [kW]
T_in = 962.9 # core inlet temp [K] (from power cycle model)
T_out = 1100 # core outlet temp [K] (from power cycle model)
T_bulk = T_in + (T_out - T_in) / 2 # bulk coolant temp. [K]
T_centerline = 1473.15 # centerline fuel temperature [K]
