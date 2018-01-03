from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys


def oneD_flow_modeling(radius, pitch, guess, L, c):
    test = FlowIteration(radius, pitch, c, L, guess)
    test.Iterate()
    test.calc_reactor_mass()
    data = test.__dict__
    print([(key, round(data[key],3)) for key in sorted(data.keys())])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("guess", type=int, help="guess N channels.")
    parser.add_argument("radius", type=float, help="cool. channel radius. [m]")
    parser.add_argument("pitch", type=float, help="fuel pitch[m]")
    parser.add_argument("core_z", type=float, help="core axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")

    args = parser.parse_args()
    
    if 2*args.radius >= args.pitch:
        print("Error: Fuel pitch must be greater than coolant channel diameter!")
        sys.exit()

    oneD_flow_modeling(args.radius, args.pitch, args.guess, args.core_z,\
            args.clad_t)
