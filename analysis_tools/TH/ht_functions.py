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
    r_channel = 0; PD = 0; c = 0; L = 0; Vol_fuel = 0; Vol_cool = 0; Vol_clad=0
    guess = 0; N_channels = 0
    mass = 0
    # temperature drop
    dt = T_centerline - T_bulk
    # heat transfer and non-dimensional coefficients
    h_bar = 0; Re = 0; Nu = 0; Pr = 0
    # flow parameters
    G_dot = 0; D_e = 0; v = 0; dp = 0
    # heat generation
    q_bar = 0; q_per_channel = 0; q_therm_check = 0
    
    def __init__(self, diameter, PD, c, L, guess):
        self.r_channel = diameter / 2.0
        self.pitch = self.r_channel * PD * 2
        self.c = c
        self.L = L
        self.guess = guess
        self.iterations = 0

    def check_converge(self):
        """Test for solution convergence and pressure limit.
        """
        if self.N_channels == self.guess:
            if self.dp < dp_allowed:
                return True
            else:
                self.N_channels += 1
                return False
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

    def get_q_bar(self, Q_therm):
        """Calculate heat transfer from fuel element. Convert q_bar max from
        cylindrical geometry to hexagonal geometry. Calculate number of required
        flow channels.
        """
        r_i = self.r_channel + self.c
        r_o = self.pitch
        # this approximation does not consider axial flux variation!!!!
        
        R_fuel = (r_o**2 / (4*k_fuel)) * ((r_i/r_o)**2 - 2*math.log(r_i/r_o) - 1)
        R_clad = (r_o/2)**2 *\
        (1-(r_i/r_o)**2)*(math.log(r_i/(r_i-self.c))/k_clad)
        R_conv = (r_o/2)**2 *\
        (1-(r_i/r_o)**2)*(1/(self.h_bar*(r_i - self.c)))
        
        q_trip_max = self.dt / (R_fuel + R_clad + R_conv)
        
        # consider axial flux variation
        A_fuel = math.sqrt(3)*self.pitch**2 / 2.0
        A_cool = (self.r_channel + self.c)**2 * math.pi
        A_fuel -= A_cool

        self.q_per_channel = 2 * q_trip_max * A_fuel * self.L
        self.q_bar = self.q_per_channel / (A_fuel * self.L)
        

        # combine actual fuel volume with conservative q-bar to estimate
        # generation per pin.
        self.N_channels = math.ceil(Q_therm / self.q_per_channel)
        
        self.Vol_fuel = A_fuel * self.L * self.N_channels
        self.Vol_cool = A_cool * self.L * self.N_channels

    def calc_dp(self):
        """Calculate pressure drop subchannel
        """
        # El Wakil (9-4)
        f = 0.184 / math.pow(self.Re, 0.2)
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = f * self.L * rho_cool * self.v * self.v / (2*self.D_e)
    
    def Iterate(self):
        """Perform Flow Calc Iteration
        """
        converge = False
        while converge == False:
            # perform necessary physics calculations
            self.mass_flux_channel()
            self.calc_nondim()
            self.get_h_bar()
            self.get_q_bar(Q_therm)
            self.calc_dp()
            converge = self.check_converge()
            self.guess = self.N_channels
            self.iterations += 1

    def calc_reactor_mass(self):
        """Based on results of the iteration, calculate the reactor mass.
        """
        self.mass = self.Vol_fuel * rho_fuel
    
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
