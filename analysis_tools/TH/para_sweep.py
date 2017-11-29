import matplotlib.pyplot as plt
import math
import copy
import numpy as np
import oned_flow as od
import argparse
from collections import OrderedDict

def sweep_PD(range, D):

    PDs = np.linspace(range[0], range[1], range[2])
    n_pins = []
    h_bar = []
    Re = []
    dP = []
    q_bar = []
    q_trip = []
    t_clad = []
    for PD_ratio in PDs:
        results = od.run_iter(PD_ratio, D)
        n_pins.append(results[0])
        h_bar.append(results[1])
        Re.append(results[2])
        dP.append(results[3])
        q_bar.append(results[4])
        t_clad.append(results[5])
        q_trip.append(results[6])
    _data = [n_pins, h_bar, Re, dP, q_bar, t_clad, q_trip]
    
    return _data, PDs

def run_analysis(minD, maxD): 
    Ds = np.linspace(minD, maxD, 5)
    for D in Ds:
        data, PDs = sweep_PD((1.0, 1.5, 40), D)
        label = str(D * 100) + ' [cm]'
        title = 'Number of Pins for: {0} <= D_f <= {1}'.format(maxD, minD)
        plot('n_pins', 1, PDs, data[0], title, 'N pins [-]', (0, 50), label, 'lower right', hline=False)
        title = 'dP Core for: {0} <= D_f <= {1}'.format(maxD, minD)
        plot('dp', 2, PDs, data[3], title, 'dp [Pa]', (0, 5e5), label, 'center right',True)
        title = 'Q_trip for: {0} <= D_f <= {1}'.format(maxD, minD)
        plot('q_trip', 3, PDs, data[4], title, "q''' [W/cc]", (0, 10), label, 'upper right', False)
        title = 'Max Clad T for: {0} <= D_f <= {1}'.format(maxD, minD)
        plot('maxT_clad', 4, PDs, data[5], title, "T_clad max [K]", (1360, 1450),
                label, 'lower right', False)
        title = 'Heat Transfer Coeff. for: {0} <= D_f <= {1}'.format(maxD, minD)
        plot('h_bar', 5, PDs, data[1], title, "h_bar [W/m^2 - K]", (0, 15), label, 'upper right', False)
        title = "Reynold's number for: {0} <= D_f <= {1}".format(maxD, minD)
        plot('re', 6, PDs, data[2], title, "Re [-]", (0, 200000), label, 'upper right', False)
       
       
def plot(savename, figure, PDs, data, title, ylabel, lim, label, loc, hline):
    plt.figure(figure)
    plt.plot(PDs, data, label=label)
    if hline:
        plt.axhline(483500, linewidth=4, label='Allowable Core dP')
    plt.title(title)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel(ylabel)
    plt.xlabel('P/D [-]')
    plt.ylim(lim[0], lim[1])
    plt.grid(which='major', alpha=1)
    plt.grid(which='minor', alpha=0.5)
    plt.xlim(min(PDs) - 0.1, max(PDs) + 0.1)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc=loc)
    plt.savefig(savename + '.png', dpi=500)
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-minD", "--min", type=float, default=0.0075,
                        required=False,
                        help="Provide min fuel diameter")
    parser.add_argument("-maxD", "--max", type=float, default=0.015, 
                        required=False,
                        help="Provide max fuel diameter")
    args = parser.parse_args()
    run_analysis(args.min, args.max)
