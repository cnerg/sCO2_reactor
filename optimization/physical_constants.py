"""Physical Constants Used for 1D simulation
*** All values are bulk flow values averaged axially across the core ***
"""
import math
import numpy as np

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """
    # default flow properties
    defaults = {'m_dot' : 0.75, # coolant flow [kg/s]
                'Q_therm' : 131000, # core thermal power [W]
                'T' : 1031.45, # bulk coolant temp [K]
                'P' : 1.766e7, # bulk coolant pressure [Pa]
                'dp_limit' : 483500, # pressure drop limit [Pa]
               }

    def __init__(self, **kwargs):
        """Inialize FlowProperties class and load required flow property data.

        Modified Attributes:
        -----------------------
            m_dot: (float) mass flow rate [kg/s]
            Q_therm: (float) required thermal reactor power [W]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_limit: (float) power-cycle constrained dp [Pa]
        """
        if kwargs:
            self.m_dot = kwargs['m_dot']
            self.Q_therm = kwargs['Q_therm']
            self.T = kwargs['T']
            self.P = kwargs['P']
            self.dp_limit = kwargs['dp_limit']
        else: # use default values
            for property in self.defaults:
                self.__dict__[property] = self.defaults[property]
            print("Warning, default flow properties are for testing purposes!!")

        # estimate secondary properties 
        self.secondary_properties()
        # if the input temperature is out of range of the fit, print a warning
        # message
        if self.T < self.t_limit[0] or self.T > self.t_limit[1]:
            print("Warning T outside of fit range. Consider re-calculating your\
 fit coeffs. to include this temperature!")

    def secondary_properties(self):
        """Calculate secondary properties from primary flow properties. Using
        two arrays of fit coefficients, perform a linear fit for all required
        flow properties. Finally, calculate Pr from mu, k_cool, and Cp

        Modified Attributes:
        --------------------
            Cp: (float) specific heat [kJ/kg-K]
            mu: (float) dynamic viscosity [kg/m-s]
            k_cool: (float) coolant conductivity [W/m-k]
            rho: (float) coolant density [kg/m^3]
            Pr: (float) cooland Prandtl number [-]
        """
        # temperature limit for curve fits
        self.t_limit = (900, 1200)
        
        # coefficients for linear fit f(T) = A*T + b
        # k, mu, rho, Cp
        A = np.array([7.0182e-5, 2.6652e-8, -0.080314, 0.255])
        B = np.array([0.0135, 1.32523e-5, 167.308, 977.66])

        # evaluate linear fit
        [self.k_cool, self.mu, self.rho, self.Cp] = np.add(A*self.T, B)
        # calculate Pr number
        self.Pr = self.Cp * self.mu / self.k_cool
         
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
