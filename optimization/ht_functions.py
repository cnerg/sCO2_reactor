# Other Imports
import math
import numpy as np
import operator
import sys
from scipy.optimize import minimize, minimize_scalar
# Import physical constants
import physical_constants as pc

def oned_flow_modeling(analyze_flow):
    """1D calculation.
    This function produces a valid, coolable reactor design given the following
    arguments:

    Arguments:
    ----------
        diameter: coolant channel diameter [m]
        PD: fuel pitch to cool. D ratio [-]
        L: reactor length [m]
        c: clad thickness [m]
    Returns:
    --------
        None
    """
    find_n_channels(analyze_flow)
    analyze_flow.adjust_dp()
    analyze_flow.calc_reactor_mass()
    analyze_flow.calc_aspect_ratio()

def _calc_n_channels_error(guess, flowiteration):
    """Calculate squared error between guess value and N channels for all
    three guess values.

    Arguments:
    ----------
        guess: a guess value for N_channels
        flowiteration: class containing TH methods and attributes

    Returns:
    --------
        error: difference between guess fuel channels and calculated required
        N_channels (float)
    """
    return flowiteration.compute_channels_from_guess(guess)


def find_n_channels(flow):
    """Perform error minimization. Using scipy's optimization package, call
    the _error function until the error is minimized to a set tolerance.

    Arguments:
    ----------
        flow: (class) Flow object. Contains attributes and
        methods required to perform an N_channels calculation for a single
        geometry (r, PD, L, c)
    Returns:
    --------
        none

    """
    res = minimize_scalar(_calc_n_channels_error, bounds=(1, 1e9), args=(flow),
                          method='Bounded', options={'xatol': 1e-3})


class Flow:
    """ Perform 1D Flow Analysis

    This class contains the required methods to perform a 1D coupled heat
    transfer/fluid flow problem on a CERMET Flow Channel.
    """
    ################################
    # UNIT SYSTEM: m, kg, J, W, Pa #
    ################################

    # geometric attributes
    r_channel = 0
    c = 0 # clad thickness
    pitch = 0 # fuel pitch (center to side of hex)
    L = 0 # reactor length
    Vol_fuel = 0 # fuel volume
    mass = 0 # fuel mass
    A_fuel = 0 # fuel cross-sectional area
    A_flow = 0 # flow cross-sectional area
    guess_channels = 0 # guess value to number of fuel channels
    N_channels = 0 # number of required fuel channels for given flow conditions
    
    # flow parameters
    Re = 0 # Reynold's number
    G_dot = 0 # mass flux
    D_e = 0 # hydraulic diameter
    v = 0 # flow velocity
    dp = 0 # channel pressure drop
    
    # heat transfer attributes
    dt = pc.T_centerline - pc.T_bulk # temp. drop fuel -> coolant
    h_bar = 0 # average heat transfer coefficient
    Nu = 0 # Nusselt number
    f = 0 # friction factor
    R_fuel = 0 # resistance to conduction in fuel
    R_clad = 0 # resistance to conduction in clad
    R_conv = 0 # resistance to convection to cool
    R_tot = 0 # total resistance to heat transfer
    
    # heat generation
    q_bar = 0 # axially-averaged volumetric generation
    q_per_channel = 0 # generation per fuel channel

    def __init__(self, diameter, PD, c, L):
        """Initialize the flow iteration class.

        Initialized Attributes:
        --------------------
            r_channel: radius of coolant channel [m]
            c: cladding thickness [m]
            pitch: fuel thickness (minor axis of hexagon) [m]
            L: length of core [m]
        """
        self.pd_ratio = PD
        self.r_channel = diameter / 2.0
        self.c = c
        self.pitch = (self.r_channel + self.c) * self.pd_ratio * 2
        self.L = L
        # set up geometry
        self.set_geom()
        # get equivalent annular radii for q_bar calculations
        self.r_i = self.r_channel + self.c
        self.r_o = self.pitch / math.sqrt(3)

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

        Correlations and equations used are from El-Wakil's Nuclear Heat
        Transport Textbook

        Modified Attributes:
        --------------------
            G_dot: mass flux through all flow channels [kg/m^2-s]
            v: flow velocity [m/s]
            Re: Reynolds number [-]
            f: friction factor [-]
            Nu: Nusselt number [-]
            h_bar: heat transfer coefficient [W/m^2-K]
        """
        # calculate mass flux
        self.G_dot = pc.m_dot / (self.A_flow * self.guess_channels)
        # calculate flow velocity from mass flux
        self.v = self.G_dot / pc.rho_cool
        # calculate Reynolds Number
        self.Re = pc.rho_cool * self.v * self.D_e / pc.mu
        # Darcy-Weisbach friction factor for pressure drop correlation El Wakil (9-4)
        self.f = 0.184 / math.pow(self.Re, 0.2)
        # Dittus-Boelter equation (9-22) from El-Wakil
        self.Nu = 0.023*math.pow(self.Re, 0.8)*math.pow(pc.Pr, 0.4)
        # heat transfer coefficient
        self.h_bar = self.Nu * pc.k_cool / self.D_e

    def get_q_per_channel(self):
        """Calculate achievable average volumetric generation:
        This method uses previously set geometry and calculated flow
        parameters to determine the maximum-achievable volumetric generation
        for each fuel channel. Calculates number of fuel channels required
        for desired thermal output.

        This method uses an equation from El-Wakil's Nuclear Heat Transport
        textbook.

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

        # El Wakil (6-62) calculates max q''' at axial centerline
        self.R_fuel = (self.r_o**2 / (4*pc.k_fuel)) *\
                      ((self.r_i/self.r_o)**2 - 2*math.log(self.r_i/self.r_o) - 1)

        self.R_clad = (self.r_o**2)/2 * (1-(self.r_i/self.r_o)**2) *\
                       math.log(self.r_i/(self.r_i-self.c)) / pc.k_clad
        
        self.R_conv = (self.r_o**2)/2 * (1-(self.r_i/self.r_o)**2) *\
                       1 / (self.h_bar*(self.r_i - self.c))

        self.R_tot = self.R_fuel + self.R_clad + self.R_conv

        # calculate centerline volumetric generation
        q_trip_max = self.dt / self.R_tot

        # consider axial flux variation
        self.q_bar = q_trip_max * 2 / math.pi
        self.q_per_channel = self.q_bar * self.A_fuel * self.L

        # calculate required fuel channels
        self.N_channels = pc.Q_therm / self.q_per_channel

    def calc_dp(self):
        """Calculate axial pressure drop across the reactor core.

        Modified Attributes:
        --------------------
            dp: core pressure drop [Pa]
        """
        # Darcy pressure drop (El-Wakil 9-3)
        self.dp = self.f * self.L * pc.rho_cool * \
                  self.v * self.v / (2*self.D_e)

    def adjust_dp(self):
        """Check for pressure constraint. This method calls calc_dp() to get
        the pressure drop in the current condition. It checks the dp against the
                self.r = radius
        self.PD = PD
power cycle-constrained allowable dp. If the pressure is too high, it
        adjusts N_channels to the min N_channels that satisfies the dp
        constraint.

        Modified Attributes:
        --------------------
            guess_channels: guess number of fuel channels [-]
            N_channels: number of fuel channels [-]
        """

        self.calc_dp()
        while self.dp > pc.dp_allowed:
            # set N_channels and guess_channels
            self.guess_channels = self.get_dp_constrained_Nchannels()
            self.N_channels = self.guess_channels
            self.characterize_flow()
            self.calc_dp()

    def get_dp_constrained_Nchannels(self):
        """Set the N_channels based on the allowable dP. This method
        calculates the required number of channels to meet the pressure drop
        constraint (set by the power cycle).

        Arguments:
        ----------
            self: Flow object [-]
        Returns:
        --------
            req_channels: Min N_channels required to meet dp constraint [-].
        """
        v_req = math.sqrt(2*self.D_e * pc.dp_allowed /
                          (self.f * self.L * pc.rho_cool))
        req_channels = math.ceil(
            pc.m_dot / (self.A_flow * pc.rho_cool * v_req))

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

    def compute_channels_from_guess(self, inp_guess):
        """Perform single 1D heat flow calculation. This method calls the
        required methods to perform one iteration of the calculation.
            
            Arguments:
            ----------
                inp_guess: (int) guess value for number of fuel channels
            Returns:
            --------
                error: (float) squared error between guess value for N_channels
                and the calculate N_channels
        """

        self.guess_channels = inp_guess
        self.characterize_flow()
        self.get_q_per_channel()

        return (self.guess_channels - self.N_channels)**2

    def calc_reactor_mass(self):
        """Based on results of the iteration, calculate the reactor mass.

        Modified Attributes:
        --------------------
            Vol_fuel: total fuel volume [m^3]
            mass: total fuel mass[kg]
        """
        self.Vol_fuel = self.A_fuel * self.L * self.N_channels
        self.mass = self.Vol_fuel * pc.rho_fuel


class ParametricSweep():
    """Class to store results of parametric sweeps for 1D flow channel analysis.

    """
    titles = {'mass': ("Total Fuel Mass", "m [kg]"),
              'Re': ("Reynolds Number", "Re [-]"),
              'G_dot': ("Reactor Mass Flux", "G_dot [kg/m^2-s]"),
              'N_channels': ("Number of Fuel Channels", "N Channels [-]"),
              'Nu': ("Nusselt Number", "Nu [-]"),
              'dp': ("Subchannel Pressure Drop", "dP [Pa]"),
              'h_bar': ("Heat Transfer Coefficient", "h [W / m^2 - K]"),
              'q_per_channel': ("Total Subchannel Generation", "q/channel [W]"),
              'q_bar': ("Average Volumetric Generation", "q_bar [W/m^3]"),
              'v': ("Flow Velocity", "v [m/s]"),
              'R_fuel': ("Resistance to Conduction in Fuel", "R_fuel [K/W]"),
              'R_clad': ("Resistance to Conduction in Clad", "R_clad [K/W]"),
              'R_conv': ("Resistance to Convection", "R_conv [K/W]"),
              'R_tot': ("Total Resistance to Heat Transfer", "R_tot [K/W]"),
              'AR': ("Approximate Core Aspect Ratio", "AR [-]")
              }

    def __init__(self, N):
        """Initialie ParametricSweep class.

        Initialized Attributes
        ----------------------
            N: (int) N^2 = number of grid points in the radius, PD mesh space.
            data: (ndarray) structured array containing results of the
            parametric sweep.
        """
        self.N = N
        # size of formats list
        N_cats = len(self.titles.keys()) + 2 # add 2 for r,pd
        self.data = np.zeros(N*N, dtype={'names' : list(self.titles.keys()) + 
                                                   ['r', 'pd'],
                                         'formats' : ['f8']*N_cats})

    def save_iteration(self, iteration, i, j, N):
        """ Save the data from each iteration of the parametric sweep. 
        """
        # 2D -> 1D index
        idx = i + j*N
        # store r, pd
        self.data[idx]['r'] = iteration.r_channel
        self.data[idx]['pd'] = iteration.pd_ratio
        for key in self.titles.keys():
            self.data[idx][key] = iteration.__dict__[key]

    def get_min_mass(self):
        """ After the parametric sweep is complete, find the minimum calculated
        fuel mass.
        """
        # search the results for minimum-mass configuration
        self.min_idx = list(self.data['mass']).index(min(self.data['mass']))
        
        # get data for min mass config
        self.min_mass = self.data[self.min_idx]['mass']
        
        # return the min_idx to get other results at minimum config
        return self.min_idx

    def disp_min_mass(self):
        """ Display the minimum mass configuration.
        """
        # get the optimal configuration
        self.minD = self.data[self.min_idx]['r']
        self.minPD = self.data[self.min_idx]['pd']
        # report the optimal configuration and it's corresponding fuel mass
        outstring = "1D Thermal Hydraulics Optimization Results:\n"
        outstring += "Min reactor mass (m = " + str(round(self.min_mass, 3))\
            + "[kg]) " + " occurs at r = " + str(self.minD) + "[m] & PD = "\
            + str(self.minPD) + "[-]."
        print(outstring)
