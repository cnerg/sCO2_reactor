from ht_functions import FlowIteration, StoreIteration
import matplotlib.pyplot as plt
import physical_constants as pc

def Iterate(r_channel, PD, c, L):

    oneD_calc = FlowIteration(r_channel, PD, c, L)
    oneD_calc.constrain_dp(500000)
    oneD_calc.mass_flux_channel()
    oneD_calc.get_h_bar()
    oneD_calc.get_q_bar(1)
    oneD_calc.calc_N_pins(pc.Q_therm)
    
    if guess_pins != oneD_calc.N_pins:
        print(oneD_calc.N_pins)
        N_pins = Iterate(r_channel, PD, c, L)

    return oneD_calc


test = Iterate(0.005, 2, 0.00031, 50)
test.calc_dp()

print(test.N_pins)
