""" Thermal-Hydraulic Optimization Module.
This module contains functions to perform a parametric sweep of reactor
geometric parameters. It calculates fuel mass and flow characteristics for each
set of valid geometric parameters.

Functions contained in this module:
    *oneD_flow_modeling
    *sweep_configs
"""
# Other imports
import argparse
import sys
# Import TH functions
from physical_constants import FlowProperties
from ht_functions import Flow, ParametricSweep
from plot import plot

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("r_lower", type=float, help="channel r lower lim [m]")
    parser.add_argument("r_upper", type=float, help="channel r upper lim [m]")
    parser.add_argument("pd_lower", type=float, help="PD lower lim [m]")
    parser.add_argument("pd_upper", type=float, help="PD upper lim [m]")
    parser.add_argument("z", type=float, help="axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")
    parser.add_argument("steps", type=int, help="parameter resolution")
    parser.add_argument("-plotkey", type=str,
                        help="parameter parameter to plot")
    parser.add_argument("-i", action='store_true', dest='show',
                        default=False, help="--display plot")

    args = parser.parse_args()

    if args.pd_lower <= 1:
        print("Error: Min fuel pitch must be greater than max coolant channel" +\
              "diameter! Set min PD > 1!")
        sys.exit()
    # setting flow properties
    primary_flow_data = {'T' : 1031.45, 
                         'P' : 1.766e7, 
                         'm_dot' : 0.75, 
                         'Q_therm' : 131000, 
                         'dp_limit' : 483500
                        }

    props = FlowProperties(flow_inputs=primary_flow_data)
    sweepresults = ParametricSweep(args.steps)
    sweepresults.sweep_geometric_configs((args.r_lower, args.r_upper),
                                         (args.pd_lower, args.pd_upper),
                                          args.z, args.clad_t, props)
    sweepresults.get_min_mass()
    sweepresults.disp_min_mass()

    if args.plotkey:
        plt = plot(sweepresults, args.plotkey, Flow.savedata)
        savename = args.plotkey + '.png'
        plt.savefig(savename, dpi=500)
        if args.show:
            plt.show()

if __name__ == '__main__':
    main()
