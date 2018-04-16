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

intpl = NearestNDInterpolator((data['mass'], data['enrich'], data['PD']), data['keff'])

print(intpl(1431, 0.87583, 1.53677))
