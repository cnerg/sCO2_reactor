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
        "AR", type=float, help="core aspect ratio [-]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")

    args = parser.parse_args()

    # perform calculation
    test = Flow(args.radius, args.clad_t, args.AR)
    oned_flow_modeling(test)
    # print results
    data = test.__dict__
    del data['fps']
    del data['fuel']
    
    for item in sorted(data): 
        disp = "{0}: {1:.3e}".format(item, data[item])
        print(disp)

    print("Power: {0:.3f} [kW]".format(test.Q_therm / 1000))

if __name__ == '__main__':
    main()
