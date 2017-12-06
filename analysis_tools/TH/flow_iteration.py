from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import physical_constants as pc
import numpy as np

def Iterate(r, PD, c, L, x, save_objects=[]):
    dp_check = False
    check_converge = False
    guess = 1
    while check_converge == False:
        flow_sweep_iteration = FlowIteration(r, PD, c, L, guess)
        
        # perform necessary physics calculations
        flow_sweep_iteration.mass_flux_channel()
        flow_sweep_iteration.calc_nondim()
        flow_sweep_iteration.get_h_bar()
        flow_sweep_iteration.get_q_bar(1)
        flow_sweep_iteration.calc_N_channels(pc.Q_therm)
        flow_sweep_iteration.calc_dp()
        print(flow_sweep_iteration.dp)
        # check N_channels for convergence
        check_converge = flow_sweep_iteration.check_converge()
        if check_converge == False:
            guess = flow_sweep_iteration.N_channels

    for item in save_objects:
        item.store_data(x, flow_sweep_iteration)  

    return flow_sweep_iteration

test = Iterate(0.01, 2, 0.0006, 0.5, 1)

print(test.__dict__)
