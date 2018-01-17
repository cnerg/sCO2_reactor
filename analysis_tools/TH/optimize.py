# Other imports
import numpy as np
import argparse
import sys
from scipy.optimize import minimize
# Import TH functions
from ht_functions import FlowIteration, ParametricSweep
import physical_constants as pc
        
def calc_mass(args, c): 
    D, PD, z = args
    flowdata = FlowIteration(D, PD, c, z)
    flowdata.Iterate()
    flowdata.calc_dp()
 
    # add pins if dp > allowable core pressure drop
    if flowdata.dp > pc.dp_allowed:
        Nchannels = flowdata.get_dp_constrained_Nchannels()
        print(Nchannels)
        valid = flowdata.check_new_channels(Nchannels)
        flowdata.N_channels = Nchannels 

    flowdata.calc_reactor_mass()

    return flowdata.mass, flowdata.N_channels
    

def get_opt_params():
    
    res = minimize(calc_mass, (0.01, 1.1, 1), bounds=((0.005, 0.1), (1.1, 5),
    (0.1,0.5)), args=(0.00031))
    print(res)

if __name__ == '__main__':
    get_opt_params()
