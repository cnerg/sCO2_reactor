from scipy import optimize
from scipy.interpolate import NearestNDInterpolator
import numpy as np
import math
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import neutronic_sweeps as ns

def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    data = np.genfromtxt(datafile, delimiter=',',
            names=ns.dimensions + ['keff', 'ave_E', 'mass'])
    
    return data

data = load_from_csv()

dp_cats = ['mass', 'PD', 'enrich']

dep_data = ()
for dim in dp_cats:
    dep_data += (data[dim],)

intpl = NearestNDInterpolator(dep_data, data['keff'])

diff = 0

for idx, line in enumerate(data):
    inp = data[idx][dp_cats]
    obs = intpl(tuple(inp))    
    err = (obs - data[idx]['keff'])**2
    diff += err


print(intpl(3525, 1.5856, 0.77))
