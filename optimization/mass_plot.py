"""Main script to run thermal-hydraulics calculations.

This script imports and utilizes the Flow class to perform a single 1D
flow calculation, producing a coolable reactor design for one set of geometric
parameters (provided as arguments).

The following functions are contained in this module:
    *oneD_flow_modeling
"""
import matplotlib.pyplot as plt
import numpy as np
# import required modules
import argparse
# import Flow class
from ht_functions import Flow, oned_flow_modeling
# import physical properties
import physical_properties as pp

critical_mass = {'UO2-CO2' : 51.07,
                 'UO2-H2O' : 51.07,
                 'UW-CO2'  : 167.2,
                 'UW-H2O'  : 167.2
                }

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

    for Q in np.arange(90000, 200000, 1000):
        x.append(Q)
        
        rxtr = make_rxtr(fuel, coolant, Q)
        
        y.append(rxtr)

    return (x, y)

def plot_results(results, dep):
    """Plot mass results as function of power
    """

    line_formats = {'UO2-CO2' : 'g--',
                    'UO2-H2O' : 'g-',
                    'UW-CO2'  : 'b--',
                    'UW-H2O'  : 'b-'}

    fig, ax = plt.subplots()

    for rxtr in sorted(results):
        y = [x.__dict__[dep] for x in results[rxtr][1]]
        power = results[rxtr][0]
        ax.plot(power, y, line_formats[rxtr], label=rxtr)
    
    plt.legend()
    plt.title('{0} vs. Thermal Power'.format(dep))
    plt.xlabel('Core Thermal Power [W]')
    plt.ylabel(dep)
    plt.savefig('{0}_vs_power.png'.format(dep), dpi=700) 

def stacked_area_plot(results):
    
    for rxtr in results:
        crit_mass = np.array([critical_mass[rxtr]]*len(results[rxtr][0]))
        fig, ax = plt.subplots()
        power = results[rxtr][0]
        mass = [x.mass for x in results[rxtr][1]]

        ax.stackplot(power, crit_mass, mass, labels=['critical mass', 'coolable mass'])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=5)
        plt.savefig('{0}_stacked_area.png'.format(rxtr)) 

def gen_data():

    rxtr_configs = ['UO2-CO2', 'UO2-H2O', 'UW-CO2', 'UW-H2O']
    results = {}

    for config in rxtr_configs:
        savefile = open(config + '.txt', 'w')
        fuel = config.split('-')[0]
        cool = config.split('-')[1]
        results[config] = power_dependence(fuel, cool)
        
        for idx, res in enumerate(results[config][0]):
            savefile.write('{0} {1}\n'.format(res,
                results[config][1][idx].Re)) 
        
        savefile.close()
    return results



data = gen_data()
plot_results(data, 'mass')
plot_results(data, 'gen_Q')
plot_results(data, 'Re')
plot_results(data, 'fuel_frac')
plot_results(data, 'core_r')
plot_results(data, 'h_bar')
stacked_area_plot(data)
