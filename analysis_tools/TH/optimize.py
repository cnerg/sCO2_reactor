# Other imports
import numpy as np
import argparse
import sys
from scipy.optimize import minimize
# Import TH functions
from ht_functions import FlowIteration, ParametricSweep
import physical_constants as pc
        
def oneD_flow_modeling(args, c):
    D = args[0]; PD = args[1]; L = args[2]
    test = FlowIteration(D, PD, c, L)
    test.oneD_calc()
    test.calc_reactor_mass()
    test.calc_aspect_ratio()
    data = test.__dict__
    
#    return abs(test.AR - 1)
    return test.mass

def get_opt_params():
    
    res = minimize(oneD_flow_modeling, (0.01, 1.1, 1), bounds=((0.005, 0.1), (1.1, 5),
    (0.1,1)), args=(0.00031))
    print(res)
    res_inputs = res['x']
    results = FlowIteration(res_inputs[0], res_inputs[1], 0.00031,\
            res_inputs[2])
    results.oneD_calc()
    results.calc_reactor_mass()
    results.calc_aspect_ratio()

    data = results.__dict__
    print([(key, str(data[key])) for key in sorted(data.keys())])

if __name__ == '__main__':
    get_opt_params()
