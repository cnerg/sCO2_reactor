"""Physical Constants Used for 1D simulation
*** All values are bulk flow values averaged axially across the core ***
"""
import math
import numpy as np

fuel_props = {'UW' : {
                'T_center'  : 1847.5, # centerline fuel temp (max) [K]
                'rho_W'     : 19300,  # clad density [kg/m^3]
                'rho_UN'    : 11300,  # fuel density [kg/m^3]
                'fuel_frac' : 0.6,    # volume fraction of fuel in CERMET
                'k_fuel'    : 51
                     },
              'UO2' : {
                'T_center'  : 1705.65, # centerline fuel temp (max) [K]
                'rho_fuel'  : 10970,  # clad density [kg/m^3]
                'k_fuel'    : 3.6      # fuel thermal conductivity [W/m-k]
                      }
             }

# mixed fuel density for Uranium Cermet Fuel
fuel_props['UW'].update( {'rho_fuel' : fuel_props['UW']['fuel_frac'] * 
                                       fuel_props['UW']['rho_UN'] +
                           (1 - fuel_props['UW']['fuel_frac'])*
                                fuel_props['UW']['rho_W']} )


class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """

    def __init__(self, coolant, power, m_dot, temp, pressure):
        """Inialize FlowProperties class and load required flow property data.

        Modified Attributes:
        -----------------------
            m_dot: (float) mass flow rate [kg/s]
            Q_therm: (float) required thermal reactor power [W]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_limit: (float) power-cycle constrained dp [Pa]
        """
        self.m_dot = m_dot
        self.Q_therm = power
        self.T = (temp[0] + temp[1]) / 2
        self.dp_limit = abs(pressure[1] - pressure[0])
        self.cool = coolant
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
        fit = {
        # temperature limit for curve fits
        'H2O' : {'t_limit' : (790, 1100),
               # coefficients for linear fit f(T) = A*T + b
               #                   k         mu         rho       Cp
               'A' : np.array([0.00015,    3.4757E-8, -0.03012,  0.2526]),
               'B' : np.array([0.0010235,  1.3711E-5, 61.25238,  2343.3392])
                },
        'CO2' : {'t_limit' : (790, 1100),
               # coefficients for linear fit f(T) = A*T + b
               #                   k           mu         rho         Cp
               'A' : np.array([5.97036E-5, 2.5034E-8, -0.062477, 0.14586]),
               'B' : np.array([0.0302958,  2.4244E-5, 134.47062, 1166.915214])
                }
        }
        # if the input temperature is out of range of the fit, print a warning
        # message
        if self.T < fit[self.cool]['t_limit'][0] or\
           self.T > fit[self.cool]['t_limit'][1]:
            print("Warning T outside of fit range. Consider re-calculating your\
 fit coeffs. to include this temperature!")
        
        # evaluate linear fit
        [self.k_cool, self.mu, self.rho, self.Cp] = fit[self.cool]['A']*self.T+\
                                                    fit[self.cool]['B']
        # calculate Pr number
        self.Pr = self.Cp * self.mu / self.k_cool
