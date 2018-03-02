"""Physical Constants Used for 1D simulation
*** All values are bulk flow values averaged axially across the core ***
"""
import math

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """
    # coefficients for linear fit f(T) = A*T + b
    coeffs = {'k'   : (7.0182e-5, 0.0135),
              'mu'  : (2.6652e-8, 1.32523e-5),
              'rho' : (-0.080314, 167.308),
              'Cp'  : (0.255, 977.66),
              # temperature limit for curve fits
              'T_limits' : (900, 1200)
             }
    
    W_e = 40000 # target electrical power [W]

    ###############################################################################
    #                                                                             #
    #                     Given Parameters From Power Cycle                       #
    #                                                                             #
    ###############################################################################
    m_dot = 0.75  # coolant flow [kg/s]
    Q_therm = 131000  # core thermal power [W]
    eta = W_e / Q_therm
    T_in = 962.9  # core inlet temp [K] (from power cycle model)
    T_out = 1100  # core outlet temp [K] (from power cycle model)
    T = T_in + (T_out - T_in) / 2  # bulk coolant temp. [K]
    rho = 86.96  # coolant density [kg/m^3]
    k_cool = 0.086  # coolant conductivity [W/m-k]
    Cp_cool = 1242  # coolant specific heat [J/kg-k]
    mu = 0.00004079  # coolant viscosity [kg/m-s]
    Pr = Cp_cool * mu / k_cool  # coolant Prandtl numbe[-]
    P_in = 1.79064e7  # inlet pressure [Pa]
    P_out = 1.74229e7  # outlet pressure [Pa]
    dp_allowed = P_in - P_out  # pressure drop limit [Pa]
    default = True

    def __init__(self, T, P, mass_flow, thermal_power):
        """Inialize FlowProperties class and load required flow property data.

        Modified Attributes:
        -----------------------
            m_dot: (float) mass flow rate [kg/s]
            Q_therm: (float) required thermal reactor power [W]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_allowed: (float) power-cycle constrained dp [Pa]
            mass: (float) required reactor mass
        """
        self.m_dot = mass_flow
        self.Q_therm = thermal_power
        self.T = float(T)
        self.P = float(P) * 1e3
        self.secondary_properties()
        self.dP_allowed = 0
        t_limit = self.coeffs['T_limits']
        # if the input temperature is out of range of the fit, print a warning
        # message
        if self.T < t_limit[0] or self.T > t_limit[1]:
            print("Warning T outside of fit range. Consider re-calculating your\
 fit coeffs. to include this temperature!")

    def secondary_properties(self):
        """Calculate secondary properties from primary flow properties
        
        Modified Attributes:
        --------------------
            Cp: (float) specific heat [kJ/kg-K]
            mu: (float) dynamic viscosity [kg/m-s]
            k_cool: (float) coolant conductivity [W/m-k]
            rho: (float) coolant density [kg/m^3]
            Pr: (float) cooland Prandtl number [-]
        """
        # calculate Cp, mu, k, rho, Pr    
        self.Cp = self.evaluate_fit('Cp')
        self.mu = self.evaluate_fit('mu')
        self.k_cool = self.evaluate_fit('k')
        self.rho = self.evaluate_fit('rho')
        self.Pr = self.Cp * self.mu / self.k_cool
    
    def evaluate_fit(self, prop):
        """Use a linear curve fit to estimate secondary flow property as a
        function of coolant temperature.

        Arguments:
        ----------
            prop: (str) a key to the coeffs used to access the coefficients to
            perform the linear fit.
        Returns:
        --------
            (float) the linear fit estimate to the property f(T)
        """
        # get coefficients for f(T) = A*T + B
        A = self.coeffs[prop][0]
        B = self.coeffs[prop][1]
        
        return A*self.T + B

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
#                            Literature Values                                #
#                                                                             #
###############################################################################
T_fuel_max = 1847.5  # centerline fuel temperature [K]
T_centerline = T_fuel_max
k_clad = 108.3  # clad conductivity: W @ 1473 K [W/m-K]
# conservative estimate for thermal conductivity at fuel centerline temperature.
k_fuel = fuel_cond(T_centerline)
rho_W = 19250  # clad density [kg/m^3]
rho_UN = 11300  # fuel density [kg/m^3]
fuel_frac = 0.6  # fraction of fuel in CERMET
# mixed density for CERMET fuel
rho_fuel = fuel_frac*rho_UN + (1-fuel_frac)*rho_W
