"""Main script to run thermal-hydraulics calculations.

This script imports and utilizes the FlowIteration class to perform a single 1D
flow calculation, producing a coolable reactor design for one set of geometric
parameters (provided as arguments).

The following functions are contained in this module:
    *oneD_flow_modeling
"""
# import required modules
import argparse
# import FlowIteration class
from ht_functions import FlowIteration

def oneD_flow_modeling(diameter, PD, L, c):
    """1D calculation.
    This function produces a valid, coolable reactor design given the following
    arguments:

    Arguments:
    ----------
        diameter: coolant channel diameter [m]
        PD: fuel pitch to cool. D ratio [-]
        L: reactor length [m]
        c: clad thickness [m]
    Returns:
    --------
        None
    """
    test = FlowIteration(diameter, PD, c, L)
    test.oneD_calc()
    test.check_dp()
    test.calc_reactor_mass()
    test.calc_aspect_ratio()

    # collect and print the results
    data = test.__dict__
    del data['error']
    print([(key, str(round(data[key], 3))) for key in sorted(data.keys())])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("diameter", type=float, help="coolant channel diameter [m]")
    parser.add_argument("PD", type=float, help="fuel pitch / cool. diameter [-]")
    parser.add_argument("core_z", type=float, help="core axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")

    args = parser.parse_args()
    
    if args.PD <= 1:
        print("Error: Fuel pitch must be greater than coolant channel\
        diameter!, set PD > 1")
        sys.exit()
    # perform calculation
    oneD_flow_modeling(args.diameter, args.PD, args.core_z,\
            args.clad_t)
