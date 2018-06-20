import operator
import math
import numpy as np
import glob
import sys

# module imports
import neutronic_sweeps as ns
import material_data as md

# sampling parameters from neutronic_sweeps
sampled_params = ns._dimensions
res = ['keff', 'mass', 'q_dens', 'BOL_U', 'EOL_U', 'rel_depl', 'dU']
names = sampled_params + res
types = ['f8']*len(names)

# conversion factors
g_to_kg = 1e-3
cc_to_l = 1e-3
kW_to_W = 1e3

def load_outputs(data_dir):
    """Load the MCNP output file string.

    Arguments:
    ----------
        data_dir (str): path to dir containing MCNP6 output files
    Returns:
    --------
        file_strings (list): list of lists of strings, 
    1 list for each output file.
    """
    file_strings = []
    files = glob.glob(data_dir)
    for output in files:
        # load file and read the lines
        fp = open(output, 'r')
        file_strings.append(fp.readlines())
    
    return file_strings

def parse_header_string(lines, data):
    """Get reactor parameters from header string of MCNP6 output file
    
    Arguments:
    ----------
        lines (list): list of strings (one for each line)
        data (ndarray): container to store results
    """
    for line in lines:
        if '1-' in line:
            vals = [float(x) for x in line.split()[1].split(',')[:6]]
            for key, val in zip(ns._dimensions, vals):
                data[key] = val
            # terminate parsing when header found
            break

def parse_keff(lines, data):
    """Parse the EOL keff from the output file.
    
    Arguments:
    ----------
        lines (list): lines of the output file (to be parsed)
        data (ndarray): container to store results
    """
    res_loc = []
    
    depletion_table_offset = 8
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_loc.append(idx)

    data['keff'] = float(lines[res_loc[-1]].split()[2])

def calc_fuel_mass_q_dens(data):
    """Calculate the reactor fuel mass and coolant density.
    
    Arguments:
    ----------
        data (ndarray): container to store results
    """
    core_r, r, PD, Q = data[['core_r', 'cool_r', 'PD', 'power']]
    # clad thickness
    c = 0.0031
    l = core_r*data['AR']
    core_v = math.pi*core_r*core_r*l

    pitch = 2*r*PD
    # calculate 'volumes' for fixed length
    v_cool = (r ** 2 * math.pi)
    # clad volume fraction
    v_clad = ((r + c)**2 - r**2)*math.pi
    # fuel volume fraction
    v_cermet = (math.sqrt(3)*pitch**2 / 2.0) - (r + c) ** 2 * math.pi 

    cell_vol = v_cool + v_clad + v_cermet
    # calculate vfracs from total cell volume
    vfrac_cermet = v_cermet / cell_vol
    
    # calculate fuel volume, power density, mass
    fuel_vol = core_v * vfrac_cermet
    data['q_dens'] = (Q / core_v) * kW_to_W / cc_to_l # W / l or MW/m^3
    data['mass'] = fuel_vol * md.rho_UN * g_to_kg
    
def parse_actinide_inventory(lines):
    """Parse actinide mass inventory.

    This function parses the MCNP6 output file and returns the actinide mass
    inventory for each step in the depletion calculation.

    Arguments:
    ----------
        lines (list): lines of the output file (to be parsed)

    Returns:
    --------
        act_inv (list): list of dicts with the actinide inventory at each time
        step.
    """
    tsteps =[]
    header_offset = 4 
    # find actinide inventory for each time step
    for idx, line in enumerate(lines):        
        if ' actinide inventory for material' in line:
            tsteps.append(idx + header_offset)

    act_inv = {}
    for i, tidx in enumerate(tsteps):
        line = lines[tidx]
        # for each time step collect actinide inventory (mass)
        while 'totals' not in line:
            data = line.split()
            ZAID = int(data[1])
            mass = float(data[2])
            # update actinide inventory and append masses 
            act_inv.setdefault(ZAID, [])
            act_inv[ZAID].append(mass)
            # next line
            tidx += 1
            line = lines[tidx]

    return act_inv

def depletion_analysis(act_inv, data, file_num):
    """Post-process the actinide inventory results.
    
    Arguments:
    ----------
        act_inv (list): list of dicts with the actinide inventory at each time
        step.
        data (ndarray): container to store results
    """
    data['BOL_U'] = act_inv[92235][0] * g_to_kg
    data['EOL_U'] = act_inv[92235][-1] * g_to_kg
    data['dU'] = data['BOL_U'] - data['EOL_U']
    data['rel_depl'] = data['dU'] / data['mass']
    
    if data['dU'] < 0:
        print("Error, non-physical depletion results, EOL uranium mass > BOL\
uranium mass. Check file {0}".format(file_num))
        sys.exit()

def save_store_data(data_dir='./results/*o'):
    """Parse every MCNP6 output for neutronics and mass/qdens results.

    Arguments:
    ----------
        data_dir (str): glob string pointing to dir with output files
    """
    # gather files
    files = glob.glob(data_dir)
    N = len(files)
    # initialize data array
    data = np.zeros(N, dtype={'names' : names, 'formats' : types})
    # unpack each MCNP6 result
    for idx, file in enumerate(files):
        print(file)
        
        fp = open(file, 'r')
        output_lines = fp.readlines()
        fp.close()
        
        # get sampled parameters from header string
        parse_header_string(output_lines, data[idx])
        # parse results from MCNP6 output file
        parse_keff(output_lines, data[idx])
        act_inv = parse_actinide_inventory(output_lines)
        # post-process the results
        calc_fuel_mass_q_dens(data[idx])
        depletion_analysis(act_inv, data[idx], file)
    
    # save data to csv file
    np.savetxt("depl_results.csv", data, delimiter=',', 
           fmt='%10.5f', header=','.join(names))
  
#if __name__ == '__main__':
#    save_store_data('/mnt/sdb/calculation_results/sa_results/*o')
#    save_store_data('./ref/*o')
