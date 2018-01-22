# Matplotlib imports
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.axis
from matplotlib import cm, rc
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter
# Other Imports
import math
import numpy as np
import operator
import sys
from scipy.optimize import minimize, minimize_scalar
# Import physical constants
from physical_constants import *
    
def _error(guess, flowiteration):
    """Calculate squared error between guess value and N channels for all
    three guess values.

    Arguments:
    ----------
        flowiteration: class containing TH methods and attributes
    
    Returns:
    --------
        error: difference between guess fuel channels and calculated required
        N_channels (float)
    """
    flowiteration.guess = guess
    flowiteration.Iterate()
    flowiteration.error = (flowiteration.guess - flowiteration.N_channels)**2
    
    return flowiteration.error

class FlowIteration:
    """ Perform 1D Flow Analysis

    This class contains the required methods to perform a 1D coupled heat
    transfer/fluid flow problem on a CERMET Flow Channel.
    """
    # geometric attributes
    r_channel = 0; c = 0; pitch = 0; L = 0; 
    Vol_fuel = 0; AR = 0; mass = 0; A_fuel = 0; A_flow = 0;
    guess = 0; N_channels = 0;
    # temperature drop
    dt = T_centerline - T_bulk
    # heat transfer and non-dimensional coefficients
    h_bar = 0; Re = 0; Nu = 0; Pr = 0; f = 0;
    R_fuel = 0; R_clad = 0; R_conv = 0; R_tot = 0;
    # flow parameters
    G_dot = 0; D_e = 0; v = 0; dp = 0
    # heat generation
    q_bar = 0; q_per_channel = 0;

    def __init__(self, diameter, PD, c, L):
        """Initialize the flow iteration class.
        
        Initialized Attributes:
        --------------------
            r_channel: radius of coolant channel [m]
            c: cladding thickness [m]
            pitch: fuel thickness (minor axis of hexagon) [m]
            L: length of core [m]
        """
        self.r_channel = diameter / 2.0
        self.c = c
        self.pitch = (self.r_channel + self.c) * PD * 2
        self.L = L

    def set_geom(self):
        """Setup the problem geometry.
        
        Modified Attributes:
        --------------------
            A_flow: flow area per fuel channel. [m^2]
            A_fuel: fuel area per fuel channel. [m^2]
            D_e: equivalent flow diameter [m]
        """
        self.A_flow = self.r_channel ** 2 * math.pi
        self.A_fuel = math.sqrt(3)*self.pitch**2 / 2.0 -\
                      (self.r_channel + self.c) ** 2 * math.pi
        self.D_e = 2.0 * self.r_channel
    
    def characterize_flow(self):
        """Calculate important non-dim and dim flow parameters. These parameters
        are required to determine generation per fuel channel.
        
        Modified Attributes:
        --------------------
            G_dot: mass flux through all flow channels [kg/m^2-s]
            v: flow velocity [m/s]
            Re: Reynolds number [-]
            Pr: Prandtl number [-]
            Nu: Nusselt number [-]
            f: friction factor [-]
            h_bar: heat transfer coefficient [W/m^2-K]
        """
        # calculate mass flux
        self.G_dot = m_dot / (self.A_flow * self.guess)
        # calculate flow velocity from mass flux
        self.v = self.G_dot / rho_cool
        # calculate Reynolds Number
        self.Re = rho_cool * self.v * self.D_e / mu
        # calculate Prandtl number
        self.Pr = Cp_cool * mu / k_cool
        # nusselt correlation is Eq (9-22) from El-Wakil Nuclear Heat Transport
        self.Nu = 0.023*math.pow(self.Re,0.8)*math.pow(self.Pr, 0.4)
        # friction factor for pressure drop correlation El Wakil (9-4)
        self.f = 0.184 / math.pow(self.Re, 0.2)
        # heat transfer coefficient 
        self.h_bar = self.Nu * k_cool / self.D_e

    def get_q_per_channel(self):
        """Calculate achievable average volumetric generation:
        This method uses previously set geometry and calculated flow
        parameters to determine the maximum-achievable volumetric generation
        for each fuel channel. Calculates number of fuel channels required
        for desired thermal output.
        
        Modified Attributes:
        --------------------
            R_fuel: resistance term (conduction in the fuel) [W/K]
            R_clad: resistance term (conduction in the clad) [W/K]
            R_conv: resistance term (convection to the fluid) [W/K]
            R_tot: total resistance to HT [W/K]
            q_per_channel: total generation in fuel channel [W]
            q_bar: axially-averaged volumetric generation in fuel [W]
            N_channels: required channels for desired Q [-]
        """
        r_i = self.r_channel + self.c
        r_o = self.pitch / math.sqrt(3)
       
        # El Wakil (9-62) calculates max q''' at axial centerline
        self.R_fuel = (r_o**2 / (4*k_fuel)) * ((r_i/r_o)**2 - 2*math.log(r_i/r_o) - 1)
        
        self.R_clad = (r_o**2)/2 * (1-(r_i/r_o)**2)*\
        math.log(r_i/(r_i-self.c)) / k_clad
        
        self.R_conv = (r_o**2)/2 * (1-(r_i/r_o)**2)*\
        1 / (self.h_bar*(r_i - self.c))
        
        self.R_tot = self.R_fuel + self.R_clad + self.R_conv
        
        # calculate centerline volumetric generation
        q_trip_max = self.dt / self.R_tot
        
        # consider axial flux variation
        self.q_per_channel = 2 * q_trip_max * self.A_fuel * self.L / math.pi
        self.q_bar = self.q_per_channel / (self.A_fuel * self.L)

    def set_N_channels(self):
        """Calculate the required number of fuel channels given the maximum heat
        generation per channel.
            
            Modified Attributes:
            --------------------
                N_channels: minium required fuel channels to meet thermal power
                requirements [-]
        """
        self.N_channels = Q_therm / self.q_per_channel

    def calc_dp(self):
        """Calculate axial pressure drop across the reactor core.
        
        Modified Attributes:
        --------------------
            dp: core pressure drop [Pa]
        """
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = self.f * self.L * rho_cool * self.v * self.v / (2*self.D_e)
        
    def check_dp(self):
        """Check for pressure constraint. This method calls calc_dp() to get
        the pressure drop in the current condition. It checks the dp against the
        power cycle-constrained allowable dp. If the pressure is too high, it
        adjusts N_channels to the min N_channels that satisfies the dp
        constraint.
        
        Modified Attributes:
        --------------------
            guess: guess number of fuel channels [-]
            N_channels: number of fuel channels [-]
        """

        self.calc_dp()
        if self.dp > dp_allowed:
            req_N = self.get_dp_constrained_Nchannels()
            # set N_channels and guess
            self.guess = req_N
            self.N_channels = req_N
            # set up new geometry and flow conditions 
            self.set_geom()
            self.characterize_flow()
            # recursive call to verifiy correct dp; accounting for f = f(Re)
            self.check_dp()
        

    def get_dp_constrained_Nchannels(self):
        """Set the N_channels based on the allowable dP. This method
        calculates the required number of channels to meet the pressure drop
        constraint (set by the power cycle).
        
        Arguments:
        ----------
            self: FlowIteration object [-]
        Returns:
        --------
            req_channels: Min N_channels required to meet dp constraint [-].
        """
        v_req = math.sqrt(2*self.D_e * dp_allowed / (self.f * self.L * rho_cool))
        req_channels = math.ceil(m_dot / (self.A_flow * rho_cool * v_req))
        
        return req_channels
    
    def calc_aspect_ratio(self):
        """Estimate the core aspect ratio (L/D).
        
        Modified Attributes:
        --------------------
            AR: core aspect ratio [-]
        """
        total_area = (self.A_fuel + self.A_flow) * self.N_channels
        equivalent_radius = math.sqrt(total_area / math.pi) 
        self.AR = self.L / (2*equivalent_radius)
    
    def Iterate(self):
        """Perform single 1D heat flow calculation. This method calls the
        required methods to perform one iteration of the calculation.
        """
        self.set_geom()
        self.characterize_flow()
        self.get_q_per_channel()
        self.set_N_channels()
        
    def oneD_calc(self):
        """Perform Flow Calc Iteration. Using scipy's optimization package, call
        the _error function until the problem is solved to a set tolerance.
        """
        res = minimize_scalar(_error, bounds=(1, 1e9), args=(self),
        method='Bounded', options={'xatol':1e-5})
        # round up to nearest N_channel
        self.guess = math.ceil(self.guess)
        self.N_channels = self.guess
        # use rounded N_channel to calculate flow characteristics
        self.set_geom()
        self.characterize_flow()
        self.get_q_per_channel()

    def calc_reactor_mass(self):
        """Based on results of the iteration, calculate the reactor mass.
        
        Modified Attributes:
        --------------------
            Vol_fuel: total fuel volume [m^3]
            mass: total fuel mass[kg]
        """
        self.Vol_fuel = self.A_fuel * self.L * self.N_channels
        self.mass = self.Vol_fuel * rho_fuel

class ParametricSweep():
    """Class to store results of parametric sweeps for 1D flow channel analysis.
    
    """
    titles = {'mass' : ("Total Fuel Mass", "m [kg]"),
              'Re' : ("Reynolds Number", "Re [-]"),
              'G_dot' : ("Reactor Mass Flux", "G_dot [kg/m^2-s]"),
              'N_channels' : ("Number of Fuel Channels", "N Channels [-]"),
              'Nu' : ("Nusselt Number", "Nu [-]"),
              'dp' : ("Subchannel Pressure Drop", "dP [Pa]"),
              'h_bar' : ("Heat Transfer Coefficient", "h [W / m^2 - K]"),
              'q_per_channel' : ("Total Subchannel Generation", "q/channel [W]"),
              'q_bar' : ("Average Volumetric Generation", "q_bar [W/m^3]"),
              'v' : ("Flow Velocity", "v [m/s]"),
              'R_fuel' : ("Resistance to Conduction in Fuel", "R_fuel [K/W]"),
              'R_clad' : ("Resistance to Conduction in Clad", "R_clad [K/W]"),
              'R_conv' : ("Resistance to Convection", "R_conv [K/W]"),
              'R_tot' : ("Total Resistance to Heat Transfer", "R_tot [K/W]"),
              'AR' : ("Approximate Core Aspect Ratio", "AR [-]")
             }

    D = 0; PD = 0;
    
    # dict to save data for plotting
    data = {k: [] for k in titles.keys()}
    min_idx = None; min_jdk = None; min_mass = 0; minD = 0; min_PD = 0;
    select_AR = False
    
    def __init__(self, D, PD, N, select_AR):
        self.D = D
        self.PD = PD
        self.N = N
        self.select_AR = select_AR
        for key in self.data:
            self.data[key] = np.empty([N,N])

    def save_iteration(self, iteration, i, j):
        """ Save the data from each iteration of the parametric sweep. 
        """
        for key in self.data.keys():
            self.data[key][i][j] = iteration.__dict__[key]

    def get_min_data(self):
        """ After the parametric sweep is complete, find the minimum calculated
        fuel mass corresponding with a valid aspect ratio. Save the flowdata for that
        point, display the PD, D and mass for the optimized configuration.
        """
        
        if False not in np.isnan(self.data['mass']):
            print("The range you have selected does not satisfy pressure drop\
 requirements, consider increasing coolant flow channel diameter")
            sys.exit()

        not_AR = 1
        if self.select_AR == True:
            not_AR = np.nan

        # calculate departure from AR = 1
        aspect_err = np.subtract(self.data['AR'], np.ones([self.N, self.N]))
        aspect_err = np.abs(aspect_err)
        # select AR error no greater than 0.33 (based on NuScale AR = 1.33)
        on_off = lambda x: not_AR if x > 0.33 else 1
        vfunc = np.vectorize(on_off)

        for key in self.data:
            self.data[key] = np.multiply(vfunc(aspect_err), self.data[key])
        
        if False not in np.isnan(self.data['AR']):
            print("The range you have selected does not satisfy aspect ratio\
 requirements, consider a new geometry range")
            sys.exit()
        
        # search the results for minimum-mass configuration
        mindices_mass = np.unravel_index(
                   np.nanargmin(self.data['mass']), self.data['mass'].shape)
        # get data for optimal AR
        self.min_idx, self.min_jdx = mindices_mass
        self.min_mass = self.data['mass'][self.min_idx][self.min_jdx]

        # save the optimal configuration
        self.minD = self.D[self.min_idx][self.min_jdx]
        self.minPD = self.PD[self.min_idx][self.min_jdx]
        
        # report the optimal configuration and it's corresponding fuel mass 
        outstring =  "1D Thermal Hydraulics Optimization Results:\n"
        outstring += "Reactor minimum mass (m = " + str(round(self.min_mass,3))\
        + "[kg]) " + " occurs at D = " + str(self.minD) +  "[m] & PD = "\
        + str(self.minPD) + "[-]."
        print(outstring)

    def plot(self, D, PD, key):
        """Produce surface plot of the flow results as function of PD and coolant
        channel diameter.
        """
        # get parametric sweep data
        M = self.data[key]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(D, PD, M, 
                               cmap=cm.viridis, linewidth=0,
                               vmin=0, vmax=np.nanmax(M),
                               antialiased=False)

        # set x/y axis labels, ticks
        ax.set_xlabel("Coolant Channel Diameter [m]", fontsize=7)
        plt.xticks(rotation=25, fontsize=6)
        ax.set_ylabel("Fuel Pitch to Coolant D Ratio [-]", fontsize=7)
        plt.yticks(rotation=25, fontsize=6)
        ax.set_zlabel(self.titles[key][1], fontsize=7)
         
        # Customize the z axis.
        ax.set_zlim(np.nanmin(M),np.nanmax(M))
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        
        # edit z tick labels
        for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(6)
        niceMathTextForm = ScalarFormatter(useMathText=True)
        ax.w_zaxis.set_major_formatter(niceMathTextForm)
        ax.ticklabel_format(axis="z", style="sci", scilimits=(0,0))
        plt.title(self.titles[key][0])
        
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5, format='%.0e')
        
        return plt
