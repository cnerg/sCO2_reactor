from scipy.optimize import curve_fit
import neutronic_sweeps as ns
import itertools
import numpy as np
import math
import matplotlib.pyplot as plt
from parse_outputs import filter_data, plot_results

def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    data = np.genfromtxt(datafile, delimiter=',',
            names=ns.dimensions + ['keff', 'ave_E', 'mass'])
    
    return data

data = load_from_csv()
#data = filter_data([('keff', 'great', 1.0)], data)

X = np.log(data['mass'])
Y = np.log(data['enrich'])
Z = np.log(data['keff'])

def log_func(data, a, b, c):
    x = a + b*data[0] + c*data[1]
    return x 

def plot_results(ind, dep, colorplot):
    """Generate Plots
    """
    # plot
    fig = plt.figure()
    plt.scatter(ind, dep, c=colorplot, s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(colorplot))))
    plt.colorbar(label='PD')
    # titles and labels
    plt.title("Error vs. Enrich (k > 1)")
    plt.savefig('figure.eps', dpi=1500, format='eps')
    
    plt.show()

    return plt

guess = (1, 1, 1)
popt, pcov = curve_fit(log_func, [X, Y], Z, p0=guess)
print(popt)
err = []

for m, e, k in zip(X, Y, data['keff']):

    fitval = math.exp(popt[0] + m*popt[1] + e*popt[2])
    err.append(abs((fitval - k) / k))

plot_results(data['mass'], err, data['PD'])
