"""Physical Constants Used for 1D simulation
*** All values are bulk flow values averaged axially across the core ***
"""
import math
import numpy as np

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
const = {'T_center' : 1847.5, # centerline fuel temp (max) [K]
         'k_clad' : 108.3,  # clad conductivity: W @ 1473 K [W/m-K]
         'rho_W' : 19250,  # clad density [kg/m^3]
         'rho_UN' : 11300,  # fuel density [kg/m^3]
         'fuel_frac' : 0.6,  # volume fraction of fuel in CERMET
         }
# mixed fuel density
const.update( {'rho_fuel' : const['fuel_frac'] * const['rho_UN'] +
                           (1 - const['fuel_frac'])*const['rho_W']} )
# conservative estimate (@centerline T) for fuel thermal conductivity
const.update( {'k_fuel' : fuel_cond(const['T_center'])} )

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """

    def __init__(self, flow_inputs=None):
        """Inialize FlowProperties class and load required flow property data.

        Modified Attributes:
        -----------------------
            m_dot: (float) mass flow rate [kg/s]
            Q_therm: (float) required thermal reactor power [W]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_limit: (float) power-cycle constrained dp [Pa]
        """
        # default flow properties
        primary_properties = {'m_dot' : 0.75, # coolant flow [kg/s]
                              'Q_therm' : 131000, # core thermal power [W]
                              'T' : 1031.45, # bulk coolant temp [K]
                              'P' : 1.766e7, # bulk coolant pressure [Pa]
                              'dp_limit' : 483500, # pressure drop limit [Pa]
                             }
        # load optional custom flow properties
        if flow_inputs:
            primary_properties = flow_inputs
        # store the flow properites
        self.m_dot = primary_properties['m_dot']
        self.Q_therm = primary_properties['Q_therm']
        self.T = primary_properties['T']
        self.P = primary_properties['P']
        self.dp_limit = primary_properties['dp_limit']

        # estimate secondary properties 
        self.secondary_properties()

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
        fit = {'t_limit' : (900, 1200),
               # coefficients for linear fit f(T) = A*T + b
               #                   k         mu         rho       Cp
               'A' : np.array([7.0182e-5, 2.6652e-8, -0.080314, 0.255]),
               'B' : np.array([0.0135,    1.32523e-5, 167.308,  977.66])
              }
        # if the input temperature is out of range of the fit, print a warning
        # message
        if self.T < fit['t_limit'][0] or self.T > fit['t_limit'][1]:
            print("Warning T outside of fit range. Consider re-calculating your\
 fit coeffs. to include this temperature!")
        
        # evaluate linear fit
        [self.k_cool, self.mu, self.rho, self.Cp] = fit['A']*self.T + fit['B']
        # calculate Pr number
        self.Pr = self.Cp * self.mu / self.k_cool
