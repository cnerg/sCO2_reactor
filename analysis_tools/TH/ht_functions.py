import math
from physical_constants import *
import matplotlib.pyplot as plt

class FlowIteration:
    """ Perform 1D Flow Analysis

    This class contains the required methods to perform a 1D coupled heat
    transfer/fluid flow problem on a CERMET Flow Channel. Upon instantiation,
    the user provides:
        flow_radius (float): radius of the coolant flow channels [m]
        PD (float): fuel cell pitch / coolant channel diameter [-]
        c (float): channel cladding thickness [m]
        L (float): fuel length [m]
        guess (int): initial number of flow channels [-]
    """
    # geometric attributes
    r_channel = 0; PD = 0; c = 0; L = 0
    guess = 0; N_channels = 0
    dt = T_centerline - T_bulk
    # heat transfer and non-dimensional coefficients
    h_bar = 0; Re = 0; Nu = 0; Pr = 0
    # flow parameters
    G_dot = 0; D_e = 0; v = 0; dp = 0
    # heat generation
    q_bar = 0; q_per_channel = 0; q_therm_check = 0
    
    def __init__(self, flow_radius, pitch, c, L, guess):
        self.r_channel = flow_radius
        self.pitch = pitch
        self.c = c
        self.L = L
        self.guess = guess
        self.iterations = 0

    def check_converge(self):
        """Test for solution convergence
        """
        
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
        flow_area = math.pi *self.r_channel**2 * self.guess 
        flow_perim = math.pi * self.r_channel * 2.0 * self.guess
        self.G_dot = m_dot / flow_area
        self.D_e = 4.0 * flow_area / flow_perim
        self.v = self.G_dot / rho_cool

    def calc_nondim(self):
        """ Calculate Reynolds number
        """
        self.Re = rho_cool * self.v * self.D_e / mu
        self.Pr = Cp_cool * mu / k_cool
        # nusselt correlation is Eq (9-22) from El-Wakil Nuclear Heat Transport
        self.Nu = 0.023*math.pow(self.Re,0.8)*math.pow(self.Pr, 0.4)
    
    def get_h_bar(self):
        """Calculate average heat transfer coefficient.
        Arguments: 
        Returns: h_bar (float) average heat transfer coefficient.
        """
        
        self.h_bar = self.Nu * k_cool / self.D_e

    def get_q_bar(self):
        """Calculate heat transfer from fuel element.
        """
        # Assume r_o that conscribes the hexagon and thus is limiting case for
        # heat transfer.
        r_o = self.pitch / math.sqrt(3)
        r_i = self.r_channel
        # Equation for flow through cylindrical element + generation provided in
        # El Wakil Nuclear Heat Transport (Eq. 5-62)
        R_cond_fuel = (1/ (4*k_f))*((r_i/r_o)**2 - 2*math.log(r_i/r_o) -1)
        R_cond_clad = 0.5*(1-(r_i/r_o)**2)*((1/k_c)*math.log(r_i / (r_i -\
            self.c)))
        R_conv = 0.5*(1-(r_i/r_o)**2)*(1/(self.h_bar*(r_i - self.c)))

        R_tot = R_cond_fuel + R_cond_clad + R_conv

        self.q_bar = 8*self.dt*(-r_o + r_i + self.c)**2 * self.L /\
               ((R_cond_fuel + 2*R_cond_clad + 2*R_conv)*r_o*r_o -\
               2*r_i*r_i*(R_cond_clad + R_conv))
        
    def calc_N_channels(self, Q_therm):
        """Calculate required number of pins based on reactor thermal power and
        q_bar
        """
        # combine actual fuel volume with conservative q-bar to estimate
        # generation per pin.
        Vol_fuel = math.sqrt(3)*self.pitch**2 * self.L / 2
        Vol_fuel -= (self.r_channel+self.c)**2 * math.pi * self.L
        self.q_per_channel = self.q_bar * Vol_fuel
        self.N_channels = math.ceil(Q_therm / self.q_per_channel)

    def calc_dp(self):
        """Calculate pressure drop subchannel
        """
        # El Wakil (9-4)
        f = 0.184 / math.pow(self.Re, 0.2)
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = f * self.L * rho_cool * self.v * self.v / (2*self.D_e)
    
    def Q_therm_check(self):
        """Check estimated thermal power (from flow calculations) against
        desired thermal power.
        """
        self.q_therm_check = self.N_channels * self.q_per_channel
       
        return self.q_therm_check

    def Iterate(self):
        """Perform Flow Calc Iteration
        """
        while self.check_converge() == False:
            # perform necessary physics calculations
            self.mass_flux_channel()
            self.calc_nondim()
            self.get_h_bar()
            self.get_q_bar()
            self.calc_N_channels(Q_therm)
            self.calc_dp()
            self.guess = self.N_channels

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
