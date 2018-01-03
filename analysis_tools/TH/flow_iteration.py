from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import physical_constants as pc
import numpy as np
import argparse
import sys

L = 1
c = 0.00031

def oneD_flow_modeling(radius, pitch, guess):
    test = FlowIteration(radius, pitch, c, L, guess)
    test.Iterate()
    data = test.__dict__
    print([(key, round(data[key],3)) for key in sorted(data.keys())])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("guess", type=int, help="guess N channels.")
    parser.add_argument("radius", type=float, help="cool. channel radius. [m]")
    parser.add_argument("pitch", type=float, help="fuel pitch[m]")

    args = parser.parse_args()
    
    if args.radius >= args.pitch:
        print("Error: Fuel pitch must be greater than coolant channel radius!")
        sys.exit()

    oneD_flow_modeling(args.radius, args.pitch, args.guess)
