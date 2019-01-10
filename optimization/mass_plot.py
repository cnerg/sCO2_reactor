"""Main script to run thermal-hydraulics calculations.

This script imports and utilizes the Flow class to perform a single 1D
flow calculation, producing a coolable reactor design for one set of geometric
parameters (provided as arguments).

The following functions are contained in this module:
    *oneD_flow_modeling
"""
# import required modules
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.optimize import curve_fit
import pandas

# import Flow class
from reactor_mass import reactor_mass
import physical_properties as pp
from ht_functions import Flow, oned_flow_modeling

labels = {'mass'        : 'Reactor Mass [kg]',
          'm_dot'       : 'Mass Flow Rate [kg/s]',
          'dp'          : 'Pressure Drop [Pa]',
          'Re'          : 'Reynolds Number [-]',
          'h_bar'       : 'Heat Transfer Coefficient [W/m^2-K]',
          'fuel_frac'   : 'Core Fuel Fraction [-]',
          'Nu'          : 'Nusselt Number [-]',
          'core_r'      : 'Core Radius [m]',
          'A_flow'      : 'Coolant Flow Area [m]',
          'gen_Q'       : 'Core Thermal Power [kW]',
          'Q_therm'     : 'Required Core Thermal Power [kW]',
          'LD'          : 'L over D [-]',
          'radius_cond' : 'Av. Distance to Conduction [m]',
          'R_cond_fuel' : 'Resistance to Conduction in Fuel [K/W]',
          'R_conv'      : 'Resistance to Convection [K/W]',
          'R_tot'       : 'Total Resistance to Heat Transfer [K/W]',
          'R_cond_clad' : 'Resistance to Conduction in Clad [K/W]',
          'fuel_mass'   : 'Fuel Mass [kg]',
          'cool_mass'   : 'Coolant Mass [kg]',
          'clad_mass'   : 'Cladding Mass [kg]',
          'refl_mass'   : 'Reflector Mass [kg]',
          'PV_mass'     : 'Pressure Vessel Mass [kg]',
          'T'           : 'Reactor Inlet Temperature [K]'
         }

T = {'CO2' : (793.8, 900),
     'H2O' : (900, 1100)
    }
P = {'CO2' : (1.7906e7, 1.7423e7), 
     'H2O' : (4.84e7, 4.77e7)
    }
m_dot = 1.2722
power = 150000

def lin_func(xdata, coeffs):
    
    return np.add(np.multiply(coeffs[0], xdata), coeffs[1])

def power_dependence(fuel, coolant, temp=None):
    """Calculate reactor mass as function of power
    """
    x = []
    y = []

    powers = np.linspace(100e3, 200e3, 500)
    data = np.zeros(len(powers), dtype={'names' : list(labels.keys()),
                                        'formats' : ['f8']*len(labels)})
    if temp:
        T_cool = (temp, T[coolant][1])
        # get coolant properties
        flow_props = pp.FlowProperties(coolant, m_dot, T_cool, P[coolant])
    else:
        # get coolant properties
        flow_props = pp.FlowProperties(coolant, m_dot, T[coolant], P[coolant])
    
    # initialize reactor model
    rxtr = Flow(0.0025, 0.00031, 1, 1, fuel, coolant, 
                'Inconel-718', 'Carbon', flow_props)
    
    for idx, Q in enumerate(powers):
        rxtr.Q_therm = Q
        rxtr.fps.m_dot = Q / (rxtr.fps.Cp * (T[coolant][1] - T[coolant][0]))
        # perform 1D calculation
        oned_flow_modeling(rxtr)
        rxtr.Q_therm /= 1e3
        rxtr.gen_Q /= 1e3
        rxtr.T = 0
        for key in labels.keys():
            data[idx][key] = rxtr.__dict__[key]
    return data

def m_dot_dependence(fuel, coolant):
    """Calculate reactor mass as function of power
    """
    x = []
    y = []

    mdots = np.linspace(100, 200, 100)
    data = np.zeros(len(mdots), dtype={'names' : list(labels.keys()),
                                        'formats' : ['f8']*len(labels)})
    
    for idx, m in enumerate(mdots):
        # get coolant properties
        flow_props = pp.FlowProperties(coolant, m, (T[0],T[1]), P[coolant])
        # initialize reactor model
        rxtr = Flow(0.0025, 0.00031, 1, power, fuel, coolant, 
                'Inconel-718', 'Carbon', flow_props)
        # perform 1D calculation
        oned_flow_modeling(rxtr)
        rxtr.T = T[0]
        
        for key in labels.keys():
            data[idx][key] = rxtr.__dict__[key]

    return data

def temp_dependence(fuel, coolant):
    """Calculate reactor mass as function of power
    """
    x = []
    y = []

    temps = np.linspace(700, 1299, 100)
    data = np.zeros(len(temps), dtype={'names' : list(labels.keys()),
                                       'formats' : ['f8']*len(labels)})
    
    for idx, Tin in enumerate(temps):
        # get coolant properties
        flow_props = pp.FlowProperties(coolant, m_dot, (Tin,T[1]), P[coolant])
        # initialize reactor model
        rxtr = Flow(0.0025, 0.00031, 1, power, fuel, coolant, 
                'Inconel-718', 'Carbon', flow_props)
        rxtr.T = Tin
        # perform 1D calculation
        oned_flow_modeling(rxtr)
        
        for key in labels.keys():
            data[idx][key] = rxtr.__dict__[key]

    return data

def plot_results(results, ind, dep):
    """Plot mass results as function of power
    """


    line_formats = {'UO2-CO2' : 'r--',
                    'UO2-H2O' : 'r-',
                    'UW-CO2'  : 'b--',
                    'UW-H2O'  : 'b-'}

    fig, ax = plt.subplots()

    for rxtr in sorted(results):
        res = results[rxtr]
        ax.scatter(res[ind], res[dep])# line_formats[rxtr], label=rxtr)
    
    plt.legend()
    title1 = " ".join(labels[dep].split()[:-1])
    title2 = " ".join(labels[ind].split()[:-1])
    plt.title('{0} vs. {1} ({2} [kg/s])'.format(title1, title2, m_dot))
    plt.xlabel(labels[ind])
    plt.ylabel(labels[dep])
    plt.savefig('{0}_vs_{1}.png'.format(dep, ind), dpi=700) 

def plot_results_temp(results, ind, dep):
    """Plot mass results as function of power
    """

    fig, ax = plt.subplots()
    y_upper = 0

    for temp in sorted(results):
        res = results[temp]
        ax.plot(res[ind], res[dep], label=int(np.average((temp, 1100))))

        if max(res['mass']) > y_upper:
            y_upper = max(res['mass'])
    
    ax.legend(title='Cool Temp [K]')

    title1 = " ".join(labels[dep].split()[:-1])
    title2 = " ".join(labels[ind].split()[:-1])
    plt.title('{0} vs. {1}'.format(title1, title2))
    plt.xlabel(labels[ind])
    plt.ylabel(labels[dep])
    plt.ylim(100, 150)
    plt.savefig('{0}_vs_{1}.png'.format(dep, ind), dpi=700) 

def gen_data():
    """Get data for all 4 reactor configurations
    """
    rxtr_configs = ['UO2-CO2', 'UW-CO2', 'UO2-H2O', 'UW-H2O']
    power_results = {}
    m_dot_results = {}
    temp_results = {}

    for config in rxtr_configs:
        fuel = config.split('-')[0]
        cool = config.split('-')[1]
                

        power_results[config] = power_dependence(fuel, cool)

    return power_results#, m_dot_results, temp_results

def gen_data_temp():
    """Get data for all 4 reactor configurations
    """
    rxtr_configs = [600, 700, 800, 900, 1000, 1100]
    power_results = {}
    m_dot_results = {}
    temp_results = {}

    for temp in rxtr_configs:

        power_results[temp] = power_dependence('UO2', 'CO2', temp)

    return power_results

def fit_power_curve(power, mass):
    """Get curve fit for power mass relationship.
    """
    popt, pcov = curve_fit(lin_func, power, mass)

    return popt, pcov

#data, mdot_data, temp_data = gen_data()
data = gen_data_temp()
#plot_results(temp_data, 'T', 'mass')
#plot_results(mdot_data, 'm_dot', 'mass')
#plot_results(mdot_data, 'm_dot', 'gen_Q')
#plot_results(data, 'gen_Q', 'dp')
plot_results(data, 'Q_therm', 'gen_Q')
plot_results_temp(data, 'gen_Q', 'mass')
#plot_results(data, 'gen_Q', 'Re')
plot_results(data, 'gen_Q', 'fuel_frac')
#plot_results(data, 'gen_Q', 'core_r')
plot_results(data, 'fuel_frac', 'core_r')
#plot_results(data, 'gen_Q', 'A_flow')
#plot_results(data, 'gen_Q', 'h_bar')
#plot_results(data, 'gen_Q', 'LD')
#plot_results(data, 'gen_Q', 'R_cond_fuel')
#plot_results(data, 'gen_Q', 'radius_cond')
#plot_results(data, 'gen_Q', 'R_conv')
#plot_results(data, 'gen_Q', 'R_tot')
#plot_results(data, 'gen_Q', 'R_cond_clad')
#plot_results(data, 'Re', 'Nu')
#stacked_area_plot(data)
#stacked_area_mass(data)
