"""Main script to run thermal-hydraulics calculations.

This script imports and utilizes the Flow class to perform a single 1D
flow calculation, producing a coolable reactor design for one set of geometric
parameters (provided as arguments).

The following functions are contained in this module:
    *oneD_flow_modeling
"""
# import required modules
import argparse
# import Flow class
from ht_functions import Flow, oned_flow_modeling

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("radius", type=float,
                        help="coolant channel radius [m]")
    parser.add_argument(
        "PD", type=float, help="fuel pitch / cool. diameter [-]")
    parser.add_argument("core_z", type=float, help="core axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")

    args = parser.parse_args()

    if args.PD <= 1:
        print("Error: Fuel pitch must be greater than coolant channel\
        diameter!, set PD > 1")
        sys.exit()
    # perform calculation
    test = Flow(args.radius, args.PD, args.clad_t, args.core_z)
    oned_flow_modeling(test)
    # print results
    data = test.__dict__
    del data['fps']
    data = {key : str(round(data[key], 3)) for key in sorted(data.keys())}
    print(data)

if __name__ == '__main__':
    main()
