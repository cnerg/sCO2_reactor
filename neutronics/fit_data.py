from scipy import optimize
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
            names=list(ns.parameters.keys()) + ['keff', 'ave_E'])
    
    return data

data = load_from_csv()

def func(data, a, b, c, d, e, f):
    """
    """
    return a*data[0]**2 + b*data[1]**2 + c*data[0]*data[1] + d*data[0] + e*data[1] * f

guess = (1,1,1,1,1,1)
C, pcov = optimize.curve_fit(func, (data['core_r'], data['PD']), data['keff'], guess)

X = data['core_r'] 
Y = data['PD']

Z = C[0]*X**2 + C[1]*Y**2 + C[2]*X*Y + C[3]*X + C[4]*Y + C[5]


print(np.mean(abs(Z - data['keff'])))
