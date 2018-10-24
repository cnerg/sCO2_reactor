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

labels = {'mass' : 'Reactor Mass [kg]',
          'Re' : 'Reynolds Number [-]',
          'h_bar' : 'Heat Transfer Coefficient [W/m^2-K]',
          'fuel_frac' : 'Core Fuel Fraction [-]',
          'Nu' : 'Nusselt Number [-]',
          'core_r' : 'Core Radius [m]',
          'A_flow' : 'Coolant Flow Area [m]',
          'gen_Q'  : 'Core Thermal Power [W]',
          'LD'     : 'L over D [-]',
          'radius_cond' : 'Av. Distance to Conduction [m]',
          'R_cond_fuel' : 'Resistance to Conduction in Fuel [K/W]',
          'R_conv' : 'Resistance to Convection [K/W]',
          'R_tot'  : 'Total Resistance to Heat Transfer [K/W]',
          'R_cond_clad' : 'Resistance to Conduction in Clad [K/W]',
          'fuel_mass' : 'Fuel Mass [kg]',
          'cool_mass' : 'Coolant Mass [kg]',
          'clad_mass' : 'Cladding Mass [kg]',
          'refl_mass' : 'Reflector Mass [kg]',
          'PV_mass'   : 'Pressure Vessel Mass [kg]',
         }

T = (900,1100)
P = {'CO2' : (1.79e7, 1.73e7), 
     'H2O' : (4.84e7, 4.77e7)
    }
m_dot = 1.15

def lin_func(xdata, coeffs):
    
    return np.add(np.multiply(coeffs[0], xdata), coeffs[1])

def power_dependence(fuel, coolant):
    """Calculate reactor mass as function of power
    """
    x = []
    y = []

    powers = np.arange(90000, 200000, 4583)
    data = np.zeros(len(powers), dtype={'names' : list(labels.keys()),
                                        'formats' : ['f8']*len(labels)})
    for idx, Q in enumerate(powers):
        rxtr = reactor_mass(fuel, coolant, Q, m_dot, T, P[coolant])
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

        ax.stackplot(res['gen_Q'], res['R_conv'], res['R_cond_fuel'],
                     res['R_cond_clad'],
                     labels=['Convection', 'Conduction Fuel', 'Conduction Clad'])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        plt.title(rxtr)
        plt.xlabel('Thermal Power [W]')
        plt.ylabel('Resistance to HT [K/W]')
        # Put a legend below current axis
        plt.legend(loc='upper right')
        plt.savefig('{0}_stacked_resistance.png'.format(rxtr)) 

def stacked_area_mass(results):
    """Plot the contributions to reactor mass for neutronic and thermal limits.
    """
    for rxtr in results:
        res = results[rxtr]
        crit_mass = [critical_mass[rxtr]]*len(results[rxtr])
        fig, ax = plt.subplots()

        ax.stackplot(res['gen_Q'], res['fuel_mass'], res['clad_mass'],
                     res['cool_mass'], res['refl_mass'], res['PV_mass'],
                     labels=['Fuel', 'Clad', 'Cool', 'Refl', 'PV'])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])
        plt.title(rxtr)
        plt.xlabel('Thermal Power [W]')
        plt.ylabel('Mass [kg]')
        # Put a legend below current axis
        plt.legend(loc='upper left')
        plt.savefig('{0}_stacked_mass.png'.format(rxtr)) 

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
                   fmt='%10.5f', header=','.join(list(labels.keys())))

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
#plot_results(data, 'gen_Q', 'A_flow')
plot_results(data, 'gen_Q', 'h_bar')
#plot_results(data, 'gen_Q', 'LD')
#plot_results(data, 'gen_Q', 'R_cond_fuel')
#plot_results(data, 'gen_Q', 'radius_cond')
#plot_results(data, 'gen_Q', 'R_conv')
#plot_results(data, 'gen_Q', 'R_tot')
#plot_results(data, 'gen_Q', 'R_cond_clad')
#plot_results(data, 'Re', 'Nu')
stacked_area_plot(data)
stacked_area_mass(data)
