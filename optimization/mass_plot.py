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
from ht_functions import Flow, oned_flow_modeling
# import physical properties
import physical_properties as pp

savevals = ['gen_Q', 'mass', 'Re', 'h_bar', 'fuel_frac', 'Nu', 'core_r', 'A_flow']
critical_mass = {'UO2-CO2' : 51.07,
                 'UO2-H2O' : 51.07,
                 'UW-CO2'  : 167.2,
                 'UW-H2O'  : 167.2
                }

def lin_func(xdata, coeffs):
    
    return np.add(np.multiply(coeffs[0], xdata), coeffs[1])

def make_rxtr(fuel, coolant, power, cool_r=0.005, clad_t=0.00031, AR=1):
    
    config = fuel + '-' + coolant
    # perform calculation
    flow_props = pp.FlowProperties(all_inp=coolant)
    rxtr = Flow(cool_r, clad_t, AR, power, config, flow_props)
    oned_flow_modeling(rxtr)
    
    return rxtr

def power_dependence(fuel, coolant):
    """Calculate reactor mass as function of power
    """
    x = []
    y = []

    powers = np.arange(90000, 200000, 400)
    

    data = np.zeros(len(powers), dtype={'names' : savevals,
                                        'formats' : ['f8']*len(savevals)})
                                    
    for idx, Q in enumerate(powers):
        rxtr = make_rxtr(fuel, coolant, Q)
        for key in savevals:
            data[idx][key] = rxtr.__dict__[key]

    return data

def plot_results(results, ind, dep):
    """Plot mass results as function of power
    """

    labels = {'mass' : 'Reactor Mass [kg]',
              'Re' : 'Reynolds Number [-]',
              'h_bar' : 'Heat Transfer Coefficient [W/m^2-K]',
              'fuel_frac' : 'Core Fuel Fraction [-]',
              'Nu' : 'Nusselt Number [-]',
              'core_r' : 'Core Radius [m]',
              'A_flow' : 'Coolant Flow Area [m]',
              'gen_Q'  : 'Core Thermal Power [W]'
             }

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
    plt.title('{0} vs. {1}'.format(title1, title2))
    plt.xlabel(labels[ind])
    plt.ylabel(labels[dep])
    plt.savefig('{0}_vs_{1}.png'.format(dep, ind), dpi=700) 

def stacked_area_plot(results):
    """Plot the contributions to reactor mass for neutronic and thermal limits.
    """
    for rxtr in results:
        res = results[rxtr]
        crit_mass = [critical_mass[rxtr]]*len(results[rxtr])
        fig, ax = plt.subplots()

        ax.stackplot(res['gen_Q'], crit_mass, res['mass'], labels=['critical mass', 'coolable mass'])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        plt.title(rxtr)
        plt.xlabel('Thermal Power [W]')
        plt.ylabel('Mass [kg]')
        # Put a legend below current axis
        plt.legend(loc='upper left')
        plt.savefig('{0}_stacked_area.png'.format(rxtr)) 

def gen_data():
    """Get data for all 4 reactor configurations
    """
    rxtr_configs = ['UO2-CO2', 'UO2-H2O', 'UW-CO2', 'UW-H2O']
    results = {}

    for config in rxtr_configs:
        fuel = config.split('-')[0]
        cool = config.split('-')[1]
        
        results[config] = power_dependence(fuel, cool)
        np.savetxt(config + '.csv', results[config], delimiter=',',
                   fmt='%10.5f', header=','.join(savevals))

    return results

def fit_power_curve(power, mass):
    """Get curve fit for power mass relationship.
    """
    popt, pcov = curve_fit(lin_func, power, mass)

    return popt, pcov

data = gen_data()
plot_results(data, 'gen_Q', 'mass')
plot_results(data, 'gen_Q', 'Re')
plot_results(data, 'gen_Q', 'fuel_frac')
plot_results(data, 'gen_Q', 'core_r')
plot_results(data, 'gen_Q', 'A_flow')
plot_results(data, 'gen_Q', 'h_bar')
plot_results(data, 'Re', 'Nu')
#stacked_area_plot(data)
