import matplotlib.pyplot as plt
import math
import copy
import numpy as np
import oned_flow as od

def sweep_PD(range, D):

    PDs = np.linspace(range[0], range[1], range[2])
    dp_data = []
    for PD_ratio in PDs:
        dp_data.append(od.run_iter(PD_ratio, D)[3])
    return dp_data, PDs

def run_analysis(): 
    Ds = np.linspace(0.0075, 0.015, 5)
    for D in Ds:
        dp_data, PDs = sweep_PD((1.1, 1.5, 25), D)
        label = str(D * 100) + ' [cm]'
        plt.plot(PDs, dp_data, label=label)
    plt.axhline(483500, linewidth=4, label='Allowable Core dP')
    plt.title('0.75 <= D_f <= 1.5')
    plt.xlabel('P/D [-]')
    plt.ylabel('dP [Pa]')
    plt.legend()
    plt.savefig('0.75_1.5.png')
    plt.show()

if __name__ == '__main__':
    run_analysis()
