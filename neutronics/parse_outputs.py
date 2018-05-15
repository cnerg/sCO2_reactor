import operator
import math
import numpy as np
import glob

# module imports
import neutronic_sweeps as ns
import material_data as md

# sampling parameters from neutronic_sweeps
sampled_params = ns._dimensions
res = ['keff', 'ave_E', 'mass', 'q_dens', 'BOL_U', 'EOL_U', 'rel_depl', 'dU']
names = sampled_params + res
types = ['f8']*len(names)

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

    Returns:
    --------
        
    """
    for line in lines:
        if '1-' in line:
            parms = line.split()[1]            
            data['AR'] = float(parms.split(',')[0])
            data['core_r'] = float(parms.split(',')[1])
            data['cool_r'] = float(parms.split(',')[2])
            data['PD'] = float(parms.split(',')[3])
            data['power'] = float(parms.split(',')[4])
            data['enrich'] = float(parms.split(',')[5])
            
            # terminate parsing when header found
            break

def parse_keff(lines, data):
    """Parse the keff data from the output file.
    
    Arguments:
    ----------
        lines (list): lines of the output file (to be parsed)
        data (ndarray): container to store results
    """
    keff = []
    err = []
    days = []
    BU = []
    
    res_loc = []
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_loc.append(idx)
        if 'print table 210' in line:
            burndx = idx + 8

    # skip the predictor calcs 
    save_res = res_loc[0::2]

    for line_num in save_res:
        keff.append(float(lines[line_num].split()[2]))
        err.append(float(lines[line_num].split()[3]))

    for burndata in lines[burndx:]:
        if burndata == '\n':
            break
        BU.append(float(burndata.split()[8]))
        days.append(float(burndata.split()[2]))
    
    data['keff'] = keff[-1]


def parse_etal(lines, data):
    """Parse energy tallies from the output file.
    Arguments:
    ----------
        tally (str): tally number
        lines (list): lines of the output file (to be parsed)
        data (ndarray): container to store results
    """
    bins = []
    vals = []
    errs = []
    tally_locations = []
    # get number of energy bins
    bupper = lines.index(' energy bins\n')
    blower = lines.index('      total bin\n')
    nbins = blower - bupper
    
    for idx, line in enumerate(lines):
        if '1tally' in line and 'nps' in line:
            tally_locations.append(idx + 11)
            
    for tally in tally_locations:
        bindata = []
        valdata = []
        errdata = []
        for idx in range(tally, tally + nbins - 1):
            bindata.append(float(lines[idx].split()[0]))
            valdata.append(float(lines[idx].split()[1]))
            errdata.append(float(lines[idx].split()[2]))
        
        bins.append(bindata)
        vals.append(valdata)
        errs.append(errdata)
    
    # calculate average neutron energy
    average = np.average(bins, weights=vals)    
    data['ave_E'] = average

def calc_fuel_mass_q_dens(data):
    """Calculate the reactor fuel mass and coolant density.
    
    Arguments:
    ----------
        data (ndarray): container to store results
    """
    g_to_kg = 0.001
    core_r, r, PD, Q = data[['core_r', 'cool_r', 'PD', 'power']]
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
    data['q_dens'] = Q / (core_v / 1000) # W / l or MW/m^3
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
    act_inv = {}
    tsteps =[]

    offset = 4

    for idx, line in enumerate(lines):        
        if ' actinide inventory for material' in line:
            tsteps.append(idx + offset)

    for i, tidx in enumerate(tsteps):
        step = 0
        line = lines[tidx + step]

        while 'totals' not in line:
            data = line.split()
            ZAID = int(data[1])
            mass = float(data[2])

            if ZAID in act_inv.keys():
                act_inv[ZAID].append(mass)
            else:
                act_inv.update({ZAID : [0]*(i) + [mass]})
            
            step += 1
            line = lines[tidx + step]

    return act_inv

def depletion_analysis(act_inv, data):
    """Post-process the actinide inventory results.
    
    Arguments:
    ----------
        act_inv (list): list of dicts with the actinide inventory at each time
        step.
        data (ndarray): container to store results
    """
    BOL = act_inv[92235][0] / 1000.0
    EOL = act_inv[92235][-1] / 1000.0
    
    data['BOL_U'] = BOL
    data['EOL_U'] = EOL

    delta_rel = abs(EOL - BOL) / data['mass']

    data['rel_depl'] = delta_rel
    data['dU'] = abs(EOL - BOL)

def save_store_data(data_dir='./results/*o'):
    """Parse every MCNP6 output for neutronics and mass/qdens results.

    Arguments:
    ----------
        data_dir (str): path to directory with MCNP6 output files
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
        string = fp.readlines()
        fp.close()
        
        # get sampled parameters from header string
        parse_header_string(string, data[idx])
        # parse results from MCNP6 output file
        parse_keff(string, data[idx])
        parse_etal(string, data[idx])
        act_inv = parse_actinide_inventory(string)
        # post-process the results
        calc_fuel_mass_q_dens(data[idx])
        depletion_analysis(act_inv, data[idx])
    
    # save data to csv file
    np.savetxt("depl_results.csv", data, delimiter=',', 
           fmt='%10.5f', header=','.join(names))
  
if __name__ == '__main__':
    save_store_data('/mnt/sdb/calculation_results/sa_results/*o')
