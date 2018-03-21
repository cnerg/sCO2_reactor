import matplotlib.pyplot as plt
import numpy as np
import glob
import neutronic_sweeps as ns

def load_outputs(data_dir):
    """Load the MCNP output file
    """
    file_strings = []
    files = glob.glob(data_dir)
    for output in files:
        # load file and read the lines
        fp = open(output, 'r')
        file_strings.append(fp.readlines())
    
    return file_strings

def parse_keff(lines):
    """Parse the keff data from the output file.
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

    return (days, BU, keff, err)

def parse_etal(tally, lines):
    """Parse energy tallies from the output file.
    """
    bins = []
    vals = []
    errs = []
    tally_locations = []
    # get number of energy bins
    bupper = lines.index(' energy bins\n')
    blower = lines.index('      total bin\n')
    nbins = blower - bupper
    
    tally_num = '{0}tally'.format(tally)

    for idx, line in enumerate(lines):
        if tally_num in line and 'nps' in line:
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

    average = np.average(bins, weights=vals)

    return (bins, vals, errs, average)

def parse_header_string(string):
    """Get important parameters from MCNP6 input header string.
    """
    for line in string:
        if '1-' in line:
            data = line.split()[1]
            AR = float(data.split(',')[0])
            PD = float(data.split(',')[1])
            cool_r = float(data.split(',')[2])
            core_r = float(data.split(',')[3])
            enrich = float(data.split(',')[4])
            power = float(data.split(',')[7])
            
            break

    return [AR, PD, cool_r, core_r, enrich, power]

def save_store_data(data_dir='./data/*'):
    """
    """
    outstrings = load_outputs(data_dir)
    names = list(ns.parameters.keys()) + ['keff', 'ave_E']
    types = ['f8']*len(names)
    N = len(outstrings)
    data = np.zeros(N, dtype={'names' : names, 'formats' : types})

    for idx, string in enumerate(outstrings):
        params = parse_header_string(string)
        data[idx]['AR'] = round(params[0], 5)
        data[idx]['PD'] = round(params[1], 5)
        data[idx]['cool_r'] = round(params[2], 5)
        data[idx]['core_r'] = round(params[3], 5)
        data[idx]['enrich'] = round(params[4], 5)
        data[idx]['power'] = round(params[5], 5)
        data[idx]['keff'] = parse_keff(string)[2][-1]
        data[idx]['ave_E'] = parse_etal('1', string)[-1]

    np.savetxt("depl_results.csv", data, delimiter=',', fmt='%10.5f',
               header=','.join(data.dtype.names))

def plot_results(data, ind, dep, colorplot=None):
    """Generate Plots
    """
    label_strings = {'AR' : 'Core Aspect Ratio[-]',
                     'PD' : 'Fuel Pitch to Coolant Channel Diameter',
                     'cool_r' : 'Coolant Channel [cm]',
                     'core_r' : 'Core Radius [cm]',
                     'enrich' : 'U-235 Enrichment [-]',
                     'power' : 'Core Thermal Power [kW]',
                     'keff' : 'k-eff [-]',
                     'ave_E' : 'average neutron energy [MeV]'
                    }
    # plot
    fig = plt.figure()
    if colorplot:
        plt.scatter(data[ind], data[dep], c=data[colorplot],
                cmap=plt.cm.get_cmap('jet', len(set(data[colorplot]))))
        plt.colorbar(ticks=list(set(data[colorplot])), label=label_strings[colorplot])
    else:
        plt.scatter(data[ind], data[dep])
    # titles and labels
    plt.title("{0} vs. {1}".format(dep, ind))
    plt.xlabel(label_strings[ind])
    plt.ylabel(label_strings[dep])

    return plt

def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    data = np.genfromtxt(datafile, delimiter=',',
            names=list(ns.parameters.keys()) + ['keff', 'ave_E'])
    
    return data

if __name__ == '__main__':
    save_store_data()
    data = load_from_csv()
    plt = plot_results(data, 'power', 'keff')
    plt.show()
