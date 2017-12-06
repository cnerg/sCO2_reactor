from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import physical_constants as pc
import numpy as np

def Iterate(r, PD, c, L, x, save_objects=[]):
    flow_sweep_iteration = FlowIteration(r, PD, c, L, 1)
    dp_check = False
    
    while dp_check == False:
        # perform necessary physics calculations
        flow_sweep_iteration.mass_flux_channel()
        flow_sweep_iteration.calc_nondim()
        flow_sweep_iteration.get_h_bar()
        flow_sweep_iteration.get_q_bar()
        flow_sweep_iteration.calc_N_pins(161)
        flow_sweep_iteration.calc_dp()
        # check pressure constraints
        if flow_sweep_iteration.dp < 400000:
            dp_check = True
        else:
            flow_sweep_iteration.guess += 1

    for item in save_objects:
        item.store_data(x, flow_sweep_iteration)  

pressure = StoreIteration("Pressure", 'dp', ('radius','dp'))
N_pins   = StoreIteration("npins", 'N_pins', ('radius','N_pins'))
h_bar    = StoreIteration("hbar", 'h_bar', ('radius','h_bar'))
save = [pressure, N_pins, h_bar]
radii = np.arange(0.001, 0.003, 0.000001)

for radius in radii:
    Iterate(radius, 2, 0.00031, 0.5, radius, save)

for item in save:
    item.plot(show=True)
