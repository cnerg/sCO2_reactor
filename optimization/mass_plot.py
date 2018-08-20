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

    for Q in np.arange(90000, 200000, 4400):
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
    plt.xlabel('Core Thermal Power [W]')
    plt.ylabel(dep)
    plt.show()



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
                results[config][1][idx].mass)) 
        
        savefile.close()
    return results

data = gen_data()
plot_results(data, 'mass')
