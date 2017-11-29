import math

m_dot = 0.75
k_c = 21e-3 # clad conductivity [kW/m-k]
k_f = 21e-3 # fuel conductivity [kW/m-k]
k_cool = 0.07527e-3 # coolant conductivity [kW/m-k], choose min of inlet/outlet k
nu = 0.00004305 # coolant viscosity [kg/m-s], chose min of inlet/outlet
Cp_cool = 1.274 # coolant specific heat [kJ/kg-k], choose min of in/out
rho_cool = 85.99 # coolant density [kg/m^3], choose min of in/out

class FlowIteration:
    
    def __init__(self, flow_radius):
        self.r_channel = flow_radius
    
    m_dot = 0.75

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

    def mass_flux_channel(self, N_elements):
        """Calculate coolant mass flux through 1 CERMET flow channel.
        Arguments: flow_radius (float) flow channel radius
                   N_elements (int) number of fuel elements
        Returns: G_dot (float): coolant mass flux
        """
        flow_area = math.pi * self.r_channel * self.r_channel 
        flow_perim = math.pi * self.r_channel * 2.0 * N_elements
        G_dot = m_dot / flow_area
        D_e = 4.0 * flow_area / flow_perim
        v = G_dot / rho_cool

        return G_dot, D_e, v
