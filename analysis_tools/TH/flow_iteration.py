from ht_functions import FlowIteration
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys


def oneD_flow_modeling(diameter, PD, guess, L, c):
    test = FlowIteration(diameter, PD, c, L, guess)
    test.Iterate()
    test.calc_reactor_mass()
    data = test.__dict__
    print([(key, round(data[key],3)) for key in sorted(data.keys())])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("guess", type=int, help="guess N channels.")
    parser.add_argument("diameter", type=float, help="coolant channel diameter [m]")
    parser.add_argument("PD", type=float, help="fuel pitch / cool. diameter [-]")
    parser.add_argument("core_z", type=float, help="core axial height [m]")
    parser.add_argument("clad_t", type=float, help="cladding thickness [m]")

    args = parser.parse_args()
    
    if args.PD <= 1:
        print("Error: Fuel pitch must be greater than coolant channel\
        diameter!, set PD > 1")
        sys.exit()

    oneD_flow_modeling(args.diameter, args.PD, args.guess, args.core_z,\
            args.clad_t)
