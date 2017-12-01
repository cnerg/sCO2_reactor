import math
from physical_constants import *

class FlowIteration:
    
    def __init__(self, flow_radius, PD, c, L):
        self.r_channel = flow_radius
        self.PD = PD
        self.c = c
        self.L = L 
        
        # constants used internally
        self.dt = T_centerline - T_bulk
        self.h_bar = 1 # 1 to avoid div by zero error
        self.Re = 0
        self.G_dot = 0
        self.D_e = 0
        self.v = 0
        self.q_bar = 0
        self.N_pins = 1
        self.dp = 0

    def get_h_bar(self):
        """Calculate average heat transfer coefficient.
        Arguments: 
        Returns: h_bar (float) average heat transfer coefficient.
        """
        
        self.Re = self.D_e * self.G_dot/ nu
        Pr = Cp_cool * nu / k_cool
        
        # equation is combination of 9-30, 9-31b in El-Wakil's Nuclear Heat
        # Transport. Valid for 1.1 <= PD <= 1.3
        Nu = (0.026 * self.PD - 0.006) * math.pow(self.Re, 0.8) * math.pow(Pr, 0.333)
        h = Nu * k_cool / self.D_e
        
        self.h_bar = h
        
    def mass_flux_channel(self):
        """Calculate coolant mass flux through 1 CERMET flow channel.
        Arguments: flow_radius (float) flow channel radius
                   self.N_pins (int) number of fuel elements
        Returns: G_dot (float): coolant mass flux
        """
        flow_area = math.pi * self.r_channel * self.r_channel * self.N_pins 
        flow_perim = math.pi * self.r_channel * 2.0 * self.N_pins
        G_dot = m_dot / flow_area
        D_e = 4.0 * flow_area / flow_perim
        v = G_dot / rho_cool
        
        self.G_dot = G_dot
        self.D_e = D_e
        self.v = v
    
    def get_q_bar(self, h_bar_switch=0):
        """Calculate heat transfer from fuel element.
        """
        r_o = self.PD*self.r_channel
        R_tot = (1/k_f) + (self.PD*self.PD - 2*math.log(self.PD) - 1) +\
                (self.PD*self.PD - 1)*((1/k_c) * math.log( (r_o + self.c) / r_o  +\
                h_bar_switch / (self.h_bar*(r_o + self.c) )))
        
        q_bar = 4 * self.dt * (r_o - self.r_channel)*(r_o - self.r_channel) * self.L /\
                (r_o*r_o*R_tot)
        self.q_bar = q_bar


    def constrain_dp(self, allowed_dp):
        """Calculate pressure drop subchannel
        Arguments: D_e (float) [m] hydraulic diameter
                   v   (float) [m/s] bulk flow velocity
        Returns:
                   dp (float) [Pa] core pressure drop
        """
        # friction factor correlation (Todreas & Kazimi eq. 9109b)
        constant = 0.1458 + 0.03632 * (self.PD - 1) - 0.03333 * (self.PD -1) ** 2
        f_s = constant / math.pow(self.Re, 0.18)

        v = f_s * self.L * rho_cool ** 2 / (self.D_e * 2 * allowed_dp)
        
        req_G_dot = v * rho_cool
        req_flow_area = m_dot / req_G_dot 
        channel_area = math.pi * self.r_channel * self.r_channel
        self.N_pins = req_flow_area / channel_area
        
    def calc_N_pins(self, Q_therm):
        """Calculate required number of pins based on reactor thermal power and
        q_bar
        """

        self.N_pins = math.ceil(Q_therm / self.q_bar)

class StoreIteration:
    """Options for Datatype are:
        h_bar
        Re
        G_dot
        D_e
        v
        q_bar
        N_pins
        dp
    """

    def __init__(self, plotname, datatype, labels):
        # containers to store iteration data
        self.data = {'name': plotname,
                     'labels' : labels,
                     'x' : [],
                     'y' : []
                    }

        self.datatype = datatype

    def store_data(self, x, iteration):
        self.data['x'].append(x)
        self.data['y'].append(iteration.__dict__[self.datatype])

    def plot(self):
        plt.plot(self.data['x'], self.data['y'], label=self.dep_variable)
        plt.xlabel(self.data['labels'][0])
        plt.ylabel(self.data['labels'][1])
        plt.title(self.data['name'])
        plt.savefig(self.plotname + '.png', dpi=500)
