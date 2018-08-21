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
#                'k_fuel'    : 51000000
                     },
              'UO2' : {
                'T_center'  : 1705.65, # centerline fuel temp (max) [K]
                'rho_fuel'  : 10970,  # clad density [kg/m^3]
                'k_fuel'    : 3.6      # fuel thermal conductivity [W/m-k]
#                'k_fuel'    : 51000000      # fuel thermal conductivity [W/m-k]
                      }

             }

# mixed fuel density for Uranium Cermet Fuel
fuel_props['UW'].update( {'rho_fuel' : fuel_props['UW']['fuel_frac'] * 
                                       fuel_props['UW']['rho_UN'] +
                           (1 - fuel_props['UW']['fuel_frac'])*
                                fuel_props['UW']['rho_W']} )


cool_props = {'CO2' : {
                  'k_cool'   : 79.082e-3,
#                  'k_cool'   : 51000000,
                  'rho'      : 233.89,
                  'mu'       : 45.905e-6,
                  'Pr'       : 0.76273,
                  'm_dot'    : 0.8, 
                  'T'        : 975,
                  'P'        : 50e6,
                  'dp_limit' : 483500
                      },

               'H2O' : {
                  'k_cool'   : 149.95e-3,
#                  'k_cool'   : 51000000,
                  'rho'      : 123.48,
                  'mu'       : 42.209e-6,
                  'Pr'       : 0.8941,
                  'm_dot'    : 0.2, 
                  'T'        : 975,
                  'P'        : 50e6,
                  'dp_limit' : 483500
                       } 
             }

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """

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
        # default flow properties
        primary_properties = {'m_dot' : 0.8, # coolant flow [kg/s]
                              'T' : 975, # bulk coolant temp [K]
                              'P' : 1.766e7, # bulk coolant pressure [Pa]
                              'dp_limit' : 483500, # pressure drop limit [Pa]
                             }
        
        if 'all_inp' in kwargs.keys():
            self.__dict__ = cool_props[kwargs['all_inp']]
            return

        # load optional custom flow properties
        if 'prime_inp' in kwargs.keys():
            primary_properties = kwargs['prime_inp']

        # store the input flow properites
        for property in primary_properties:
            self.__dict__[property] = primary_properties[property]
        
        # load secondary properties
        if 'second_inp' in kwargs.keys():
            for prop in kwargs['second_inp']:
                self.__dict__[property] = kwargs['second_inp'][prop]
        else:
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
