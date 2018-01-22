""" Thermal-Hydraulic Optimization Module.
This module contains functions to perform a parametric sweep of reactor
geometric parameters. It calculates fuel mass and flow characteristics for each
set of valid geometric parameters.

Functions contained in this module:
    *oneD_flow_modeling
    *sweep_configs
"""
# Other imports
import numpy as np
import argparse
import sys
# Import TH functions
from ht_functions import FlowIteration, ParametricSweep

def oneD_flow_modeling(analyze_flow):
    """Conduct oneD_flow_modeling.

    Arguments:
    ----------
        analyze_flow: (class) FlowIteration object. Contains attributes and
        methods required to perform an N_channels calculation for a single
        geometry (r, PD, L, c)

    Returns:
    --------
        None
    """
    analyze_flow.oneD_calc()
    analyze_flow.check_dp()
    analyze_flow.calc_reactor_mass()
    analyze_flow.calc_aspect_ratio()
        
def sweep_configs(D, PD, z, c, N, key, AR_select, save=False):
    """Perform parametric sweep through pin cell geometric space.

    """
    # calculate appropriate step sizes given range
    D_step = (D[1] - D[0]) / N
    PD_step = (PD[1] - PD[0]) / N
    # ranges for diameter and pitch/diameter ratio
    D = np.arange(D[0], D[1], D_step)
    PD = np.arange(PD[0], PD[1], PD_step)

    # create parameter mesh
    D, PD = np.meshgrid(D, PD)
    
    # initialize object to save sweep results
    sweepresults = ParametricSweep(D, PD, N, AR_select)
    # sweep through parameter space, calculate min mass
    for i in range(N):
        for j in range(N):
            flowdata = FlowIteration(D[i,j], PD[i,j], c, z)
            oneD_flow_modeling(flowdata)
            sweepresults.save_iteration(flowdata, i, j)
    
    sweepresults.get_min_data()
    plt = sweepresults.plot(D, PD, key)
    
    if save == True:
        savename = key + '.png'
        plt.savefig(savename, dpi=500)

    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("d_lower", type=float, help="channel D lower lim [m]")
    parser.add_argument("d_upper", type=float, help="channel D upper lim [m]")
    parser.add_argument("pd_lower", type=float, help="PD lower lim [m]")
    parser.add_argument("pd_upper", type=float, help="PD upper lim [m]")
    parser.add_argument("z", type=float, help="axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")
    parser.add_argument("steps", type=int, help="parameter resolution")
    parser.add_argument("plotkey", type=str, help="parameter parameter to plot")
    parser.add_argument("--AR", action='store_true', dest='AR_select',
            default=False, help="--selects data corresponding to valid AR's")

    args = parser.parse_args()
    
    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel\
diameter! Set min PD > 1!")
        sys.exit()

    sweep_configs((args.d_lower, args.d_upper),
                                   (args.pd_lower, args.pd_upper), 
                                   args.z, args.clad_t, args.steps,
                                   args.plotkey, args.AR_select, True)
