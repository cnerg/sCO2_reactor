"""Physical Constants Used for 1D simulation of sCO2 reactor
*** All values are bulk flow values averaged axially across the core ***
"""
import numpy as np
from scipy import interpolate

kJ_to_J = 1e3

fuel_props = {'UW' : {
                'T_center'  : 1847.5, # centerline fuel temp (max) [K]
                'rho_W'     : 19300,  # clad density [kg/m^3]
                'rho_UN'    : 11300,  # fuel density [kg/m^3]
                'fuel_frac' : 0.6,    # volume fraction of fuel in CERMET
                'k'         : 51
                     },
              'UO2' : {
                'T_center'  : 1705.65, # centerline fuel temp (max) [K]
                'rho_fuel'  : 10970,  # clad density [kg/m^3]
                'k'         : 3.6      # fuel thermal conductivity [W/m-k]
                      }
              }

clad_props = {'Inconel-718' :  {
                 # assume SS for roughness of Inconel Cladding
                'rough'     : 1.5e-6,
                # from ESPI Metals Website
                'k'         : 11.4,
                'rho'       : 8.19
                      }
             }

refl_props = {'Carbon'      :  {
                'rho'       : 1700
                      }
             }

# mixed fuel density for Uranium Cermet Fuel
fuel_props['UW'].update( {'rho_fuel' : fuel_props['UW']['fuel_frac'] * 
                                       fuel_props['UW']['rho_UN'] +
                           (1 - fuel_props['UW']['fuel_frac'])*
                                fuel_props['UW']['rho_W']} )


def load_data(coolant):
    """Load coolant property data from EES tables.

    Arguments:
    ----------
        coolant (str): type of coolant
    
    Returns:
    --------
        data (ndarray): flow property data
    """
    # open datafile
    filename = './data/{0}.txt'.format(coolant)
    lines = open(filename, 'r').readlines()
    
    data = np.zeros(len(lines), dtype={
                               'names' : ['T', 'P', 'Cp', 'rho', 'k', 'mu'],
                               'formats' : ['f8']*6})
    # collect table data
    for idx, line in enumerate(lines):
        res = [float(x) for x in line.split()]
        data[idx]['T']   = res[0]
        data[idx]['P']   = res[1]
        data[idx]['Cp']  = res[2]*kJ_to_J
        data[idx]['rho'] = res[5]
        data[idx]['k']   = res[3]
        data[idx]['mu']  = res[4]

    return data

def make_interps(coolant):
    """Create interpolation functions from property data.

    Arguments:
    ----------
        coolant (str): type of coolant
    
    Returns:
    --------
        interp_props (dict): library of interpolation objects for estimating
        flow properties
        t_range (tuple): temperature range for interpolation validity
        p_range (tuple): pressure range for interpolation validity
    """
    # load data from EES tables
    data = load_data(coolant)
    points = list(zip(data['T'], data['P']))
    props = ['Cp', 'rho', 'k', 'mu']
    
    # get ranges for interpolation validity check
    t_range = (min(data['T']), max(data['T']))
    p_range = (min(data['P']), max(data['P']))

    interp_props = {}
    
    # generate the interpolation objects for each property
    for p in props: 
        interp_props[p] = interpolate.LinearNDInterpolator(points, data[p])
    
    return interp_props, t_range, p_range

def bulk_average(vals):
    """Bulk average primary property

    Arguments:
    ----------
        inlet (float): inlet property value
        outlet (float): outlet property value
    Returns:
    --------
        averaged (float): bulk averaged value
    """
    inlet = vals[0]
    outlet = vals[1]

    return (inlet + outlet) / 2

class FlowProperties:
    """Class to store flow properties and calculate secondary properties from
    fundamental props.
    """

    def __init__(self, coolant, m_dot, temp, pressure):
        """Inialize FlowProperties class and load required flow property data.

        Modified Attributes:
        -----------------------
            m_dot: (float) mass flow rate [kg/s]
            T: (float) bulk coolant temperature [K]
            P: (float) bulk coolant pressure [Pa]
            dp_limit: (float) power-cycle constrained dp [Pa]
        """
        self.m_dot = m_dot
        self.cool = coolant
        # bulk-averaged properites
        self.T = bulk_average(temp)
        self.P = bulk_average(pressure)
        self.dp_limit = abs(pressure[1] - pressure[0])
        # interpolate secondary properties
        self.prop_interps, self.t_limit, self.p_limit = make_interps(self.cool)
        self.secondary_properties()
    
    def update_props(self, m_dot, temp, pressure):
        """Update the properties with a new temperature and pressure
        """
        self.m_dot = m_dot
        self.T = bulk_average(temp)
        self.P = bulk_average(pressure)
        # interpolate new secondary properties
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
        # if the input temperature is out of range of the fit, print a warning
        # message
        if self.T < self.t_limit[0] or self.T > self.t_limit[1]:
            print("Warning T outside of fit range. Consider re-calculating your\
 fitted data to include this temperature!")
        if self.P < self.p_limit[0] or self.P > self.p_limit[1]:
            print("Warning P outside of fit range. Consider re-calculating your\
 fitted data to include this pressure!")
        
        # evaluate linear interpolator for each property
        self.k_cool = self.prop_interps['k'](self.T, self.P)
        self.mu     = self.prop_interps['mu'](self.T, self.P)
        self.rho    = self.prop_interps['rho'](self.T, self.P)
        self.Cp     = self.prop_interps['Cp'](self.T, self.P)
        # calculate Pr number
        self.Pr = self.Cp * self.mu / self.k_cool
