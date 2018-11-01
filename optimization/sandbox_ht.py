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
          'v'           : 'Flow Velocity [m/s]',
          'Re'          : 'Reynolds Number [-]',
          'G_dot'       : 'Mass Flux [kg/m^2-s]',
          'h_bar'       : 'Heat Transfer Coefficient [W/m^2-K]',
          'fuel_frac'   : 'Core Fuel Fraction [-]',
          'Nu'          : 'Nusselt Number [-]',
          'core_r'      : 'Core Radius [m]',
          'A_flow'      : 'Coolant Flow Area [m]',
          'gen_Q'       : 'Core Thermal Power [W]',
          'Q_therm'     : 'Required Core Thermal Power [W]',
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
         }

T = (900,1100)
P = {'CO2' : (1.79e7, 1.73e7), 
     'H2O' : (4.84e7, 4.77e7)
    }
m_dot = 0.2
power = 200000


def frac_dependence(fuel, coolant):
    """Calculate reactor mass as function of power
    """
    x = []
    y = []

    fracs = np.linspace(0.4, 0.99, 100)
    data = np.zeros(len(fracs), dtype={'names' : list(labels.keys()),
                                        'formats' : ['f8']*len(labels)})
    # get coolant properties
    flow_props = pp.FlowProperties(coolant, m_dot, (T[0],T[1]), P[coolant])
    # initialize reactor model
    rxtr = Flow(0.005, 0.00031, 1, power, fuel, coolant, 
                'Inconel-718', 'Carbon', flow_props)
    
    for idx, frac in enumerate(fracs):
        
        # perform 1D calculation
        rxtr.compute_Q_from_guess(frac)
        rxtr.adjust_dp()
        rxtr.calc_reactor_mass()
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
        ax.plot(res[ind], res[dep], line_formats[rxtr], label=rxtr)
    
    plt.legend()
    title1 = " ".join(labels[dep].split()[:-1])
    title2 = " ".join(labels[ind].split()[:-1])
    plt.title('{0} vs. {1} ({2} [kg/s]'.format(title1, title2, m_dot))
    plt.xlabel(labels[ind])
    plt.ylabel(labels[dep])
    plt.savefig('{0}_vs_{1}.png'.format(dep, ind), dpi=700) 

def gen_data():
    """Get data for all 4 reactor configurations
    """
    rxtr_configs = ['UO2-CO2', 'UO2-H2O', 'UW-CO2', 'UW-H2O']
#    rxtr_configs = ['UW-CO2']
    power_results = {}

    for config in rxtr_configs:
        print(config)
        fuel = config.split('-')[0]
        cool = config.split('-')[1]
        
        power_results[config] = frac_dependence(fuel, cool)

    return power_results

data = gen_data()
plot_results(data, 'fuel_frac', 'gen_Q')
plot_results(data, 'fuel_frac', 'Re')
plot_results(data, 'fuel_frac', 'A_flow')
plot_results(data, 'fuel_frac', 'G_dot')
plot_results(data, 'fuel_frac', 'v')

plot_results(data, 'Re', 'Nu')
