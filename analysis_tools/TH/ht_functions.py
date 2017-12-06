import math
from physical_constants import *
import matplotlib.pyplot as plt


class FlowIteration:
    """ Perform 1D Flow Analysis

    This class contains the required methods to perform a 1D coupled heat
    transfer/fluid flow problem on a CERMET Flow Channel. Upon instantiation,
    the user provides:
        flow_radius (float): radius of the coolant flow channels
        PD (float): fuel cell pitch / coolant channel diameter
        c (float): channel cladding thickness
        guess (int): an initial guess for the required number of flow channels 

    """
    
    def __init__(self, flow_radius, PD, c, L, guess):
        self.r_channel = flow_radius
        self.PD = PD
        self.c = c
        self.L = L
        self.guess = guess
        self.dt = T_centerline - T_bulk
        self.h_bar = 1 # 1 to avoid div by zero error
        self.Re = 0
        self.Nu = 0
        self.Pr = 0
        self.G_dot = 0
        self.D_e = 0
        self.v = 0
        self.q_bar = 0
        self.N_channels = 0
        self.dp = 0

    def check_converge(self):
        if self.N_channels == self.guess:
            return True
        else:
            return False
    
    def mass_flux_channel(self):
        """Calculate coolant mass flux through 1 CERMET flow channel.
        Arguments: flow_radius (float) flow channel radius
                   self.N_channels (int) number of fuel elements
        Returns: G_dot (float): coolant mass flux
        """
        flow_area = math.pi *(self.r_channel-self.c)**2 * self.guess 
        flow_perim = math.pi * self.r_channel * 2.0 * self.guess
        self.G_dot = m_dot / flow_area
        self.D_e = 4.0 * flow_area / flow_perim
        self.v = self.G_dot / rho_cool

    def calc_nondim(self):
        """ Calculate Reynolds number
        """
        self.Re = rho_cool * self.v * self.L / mu
        self.Pr = Cp_cool * mu / k_cool
        # nusselt correlation is Eq (9-22) from El-Wakil Nuclear Heat Transport
        self.Nu = 0.023*math.pow(self.Re,0.8)*math.pow(self.Pr, 0.4)
    
    def get_h_bar(self):
        """Calculate average heat transfer coefficient.
        Arguments: 
        Returns: h_bar (float) average heat transfer coefficient.
        """
        
        self.h_bar = self.Nu * k_cool / self.D_e
        
    def get_q_bar(self, h_bar_switch=0):
        """Calculate heat transfer from fuel element.
        """
        r_o = self.PD*self.r_channel*2 / math.sqrt(3)
        R_tot = (1/k_f) + (self.PD*self.PD - 2*math.log(self.PD) - 1) +\
                (self.PD*self.PD - 1)*((1/k_c) * math.log( (r_o + self.c) / r_o  +\
                h_bar_switch / (self.h_bar*(r_o + self.c) )))
        
        self.q_bar = 4 * self.dt * (r_o - self.r_channel)*(r_o - self.r_channel) * self.L /\
                (r_o*r_o*R_tot)

    def calc_N_channels(self, Q_therm):
        """Calculate required number of pins based on reactor thermal power and
        q_bar
        """
        pitch = self.r_channel * self.PD / 2.0
        Vol_fuel = 2*math.sqrt(3)*self.r_channel*self.r_channel * self.L
        Vol_fuel -= self.r_channel*self.r_channel * math.pi * self.L
        q_per_channel = self.q_bar * Vol_fuel
        self.N_channels = math.ceil(Q_therm / q_per_channel)

    def calc_dp(self):
        """Calculate pressure drop subchannel
        """
        # El Wakil (9-4)
        f = 0.184 / math.pow(self.Re, 0.2)
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = f * self.L * rho_cool * self.v * self.v / (2*self.D_e)

class StoreIteration:
    """ Save Data For Analysis:
    
    This object is used to store the results of the FlowIterations (see class
    above). Upon instantiation, user must provide:

    plotname: (str) Name of the plot, also used to save the plot.png file.
    datatype: (str) identifier of the data the user wants to save.
    labels: (tuple(str, str) ) desired plot labels

    Options for datatypes are:
        h_bar (average heat transfer coefficient)
        Re (reynolds number for coolant flow)
        G_dot (coolant mass flux)
        D_e (hydraulic equivalent diameter)
        v (coolant velocity)
        q_bar (average heat generation per fuel cell)
        N_channels (number of coolant channels)
        dp (pressure loss through the channels)
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
        """Store the desired data (see datatype) and the provided independent
        variable.
        """
        self.data['x'].append(x)
        self.data['y'].append(iteration.__dict__[self.datatype])

    def plot(self, show=False):
        """Plot the stored data
        """
        plt.plot(self.data['x'], self.data['y'])
        plt.xlabel(self.data['labels'][0])
        plt.ylabel(self.data['labels'][1])
        plt.title(self.data['name'])
        plt.savefig(self.data['name'] + '.png', dpi=500)
        # If user provides option show=True, display plot
        if show == True:
            plt.show()
