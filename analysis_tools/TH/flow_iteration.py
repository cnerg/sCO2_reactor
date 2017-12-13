from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import physical_constants as pc
import numpy as np
import argparse

L = 0.5
c = 0.00031

def oneD_flow_modeling(radius, PD, guess):
    test = FlowIteration(radius, PD, c, L, guess)
    test.Iterate()
    data = test.__dict__
    print([(key, data[key]) for key in sorted(data.keys())])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("guess", type=int, help="guess N channels.")
    parser.add_argument("radius", type=float, help="cool. channel radius. [m]")
    parser.add_argument("PD", type=float, help="fuel pitch to coolant channel\
    diameter ratio.")

    args = parser.parse_args()

    oneD_flow_modeling(args.radius, args.PD, args.guess)
