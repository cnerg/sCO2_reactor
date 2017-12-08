from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import physical_constants as pc
import numpy as np

def Iterate(r, PD, c, L, guess):
    iteration = FlowIteration(r, PD, c, L, guess)        
    # perform necessary physics calculations
    iteration.mass_flux_channel()
    iteration.calc_nondim()
    iteration.get_h_bar()
    iteration.get_q_bar(1)
    iteration.calc_N_channels(pc.Q_therm)
    iteration.calc_dp()
    
    if iteration.check_converge() == False:
        iteration =\
        Iterate(r, PD, c, L, iteration.N_channels)

    return iteration


test = Iterate(0.002, 100, 0.00031, 0.5, 100)
print(test.__dict__)


"""


q_bar1 = []
q_bar2 = []
q_bar3 = []
pitch = []

ratios = np.arange(1,100,1)

for PD in ratios:

    test1 = Iterate(0.01, PD, 0.00031, 0.5, 100)
    test2= Iterate(0.1, PD, 0.00031, 0.5, 100)
    test3 = Iterate(0.3, PD, 0.00031, 0.5, 100)
    pitch.append(PD)
    q_bar1.append(test1.N_channels)
    q_bar2.append(test1.N_channels)
    q_bar3.append(test1.N_channels)


plt.plot(pitch, q_bar1, pitch, q_bar2, pitch, q_bar3)
plt.show()
"""
