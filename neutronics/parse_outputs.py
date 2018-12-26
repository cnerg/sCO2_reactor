import sys
import zipfile
import operator
import math
import material_data as md
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
import neutronic_sweeps as ns
import pandas
from mcnp_inputs import HomogeneousInput

names =  ['core_r', 'fuel_frac', 'ref_mult', 'keff', 
          'stdv', 'r_mass', 'c_mass', 'p_mass', 'mass']
types = ['f8']*len(names)

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

def parse_keff(lines, filename):
    """Parse the output for keff value.
    """

    res_idx = []
    for idx, line in enumerate(lines):
        if 'final result' in line:
            res_idx.append(idx)

    keff = float(lines[res_idx[-1]].split()[2])
    stdv = float(lines[res_idx[-1]].split()[3])

    return keff, stdv

def calc_mass(data):
    """
    """
    config = {'fuel' : 'UO2',
              'matr' : None,
              'cool' : 'CO2',
              'clad' : 'Inconel-718',
             }
    config['core_r'] = data['core_r']
    config['fuel_frac'] = data['fuel_frac']
    config['ref_mult'] = data['ref_mult']

    input = HomogeneousInput(config=config)
    input.homog_core()
    
    data['mass'] = input.tot_mass
    data['r_mass'] = input.refl_mass
    data['p_mass'] = input.PV_mass
    data['c_mass'] = input.core_mass 

def parse_header_string(lines, data):
    """
    """
    for line in lines:
        if '1-' in line:
            break
    r, f, m = [float(x) for x in line.split()[1].split(',')[:-1]]

    data['core_r'] = r
    data['fuel_frac'] = f
    data['ref_mult'] = m

def save_store_data(data_dir, name='crit_results.csv'):
    """
    """
    # load the results from zip file
    zip_folder = zipfile.ZipFile(data_dir)
    files = zip_folder.namelist()
    N = len(files)
    data = np.zeros(N, dtype={'names' : names, 'formats' : types})
    

    for idx, file in enumerate(files):
        print(file)
        fp = zip_folder.open(file)
        string = [x.decode('ascii') for x in fp.readlines()]
        fp.close()
        parse_header_string(string, data[idx])
        calc_mass(data[idx])
        data[idx]['keff'], data[idx]['stdv'] = parse_keff(string, file)

    zip_folder.close()
    
    np.savetxt(name, data, delimiter=',', 
           fmt='%10.5f', header=','.join(names), comments='')
  
def plot_results(data, ind, dep, plt, colorplot=None, log=None):
    """Generate Plots
    """
    label_strings = {'core_r' : 'Core Radius [cm]',
                     'enrich' : 'U-235 Enrichment [-]',
                     'keff' : 'k-eff [-]',
                     'mass' : 'reactor fuel mass [kg]',
                     'fuel_frac' : 'volume fraction fuel [-]',
                     'ref_mult' : 'Reflector Multiplier [-]',
                     'r_mass' : 'Reflector Mass [kg]'
                    }
    # plot
#    fig = plt.figure()
    if colorplot:
        colorsave = '_'+colorplot
        plt.scatter(data[ind], data[dep], c=data[colorplot], s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(data[colorplot]))))
        plt.colorbar(label=label_strings[colorplot])
    else:
        colorsave = ''
        plt.scatter(data[ind], data[dep], s=6)
    # titles and labels
    plt.title("{0} vs. {1}".format(dep, ind))
#    plt.title(r'EOL k$_{eff}$ vs. Fuel Mass')
#    plt.title("keff vs. mass for 0.2 < enrich < 0.3")
    plt.xlabel(label_strings[ind])
    plt.ylabel(label_strings[dep])
#    plt.ylim(0.6, 1.8)

    if log == 'log-log':
        plt.xscale('log')
        plt.yscale('log')
    if log == 'semilogy':
        plt.yscale('log')
    if log == 'semilogx':
        plt.xscale('log')

    savename = '{0}_vs_{1}{2}.png'.format(dep, ind, colorsave)
#    plt.savefig(savename, dpi=1000, format='png')

    return plt


def surf_plot(data):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    X = data['core_r']
    Y = data['fuel_frac']
    Z = data['ref_mult']
    k = data['keff']

    # Plot the surface.
    ax.scatter(X,Y,Z, c=k)
    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    plt.show()

def interpolate_grid(data):
    
    skip_dx = 1
    X = sorted(list(set(data['core_r'])))[0::skip_dx]
    Y = sorted(list(set(data['fuel_frac'])))[0::skip_dx]
    Z = sorted(list(set(data['ref_mult'])))[0::skip_dx]
    keff = data['keff']
    
    kk = np.zeros([len(X), len(Y), len(Z)])
    for i, r in enumerate(X):
        for j, f in enumerate(Y):
            for k, m in enumerate(Z):
                K = data[(data['core_r'] == r) &\
                         (data['fuel_frac'] == f) &\
                         (data['ref_mult'] == m)]['keff']
                stdv = data[(data['core_r'] == r) &\
                         (data['fuel_frac'] == f) &\
                         (data['ref_mult'] == m)]['stdv']
                kk[i, j, k] = K
    
    fn = RegularGridInterpolator((X, Y, Z), kk, method='linear')
    
    return fn 

def test_interpolator(data, func, config):
    """
    """
    error = []
    for idx, row in data.iterrows():
        r = row['core_r']
        f = row['fuel_frac']
        m = row['ref_mult']
        if r > func.grid[0][-1] or r < func.grid[0][0]:
            continue
        if f > func.grid[1][-1] or f < func.grid[1][0]:
            continue
        if m > func.grid[2][-1] or m < func.grid[2][0]:
            continue

        kinterp = func([row['core_r'], row['fuel_frac'], row['ref_mult']])[0]
        error.append(abs(kinterp - row['keff']) * 1e5)
    print(len(error))
    fig = plt.figure()
    plt.hist(error)
    title = ' '.join(config.split('_'))
    plt.title('Interpolation Error: ' + title)
    plt.xlabel('Interpolated Reactivity Error [pcm]')
    plt.savefig('check_interp.png')


def load_from_csv(datafile="depl_results.csv"):
    """load the results data from a csv.
    """
    data = pandas.read_csv(datafile)
    
    return data

def filter_data(filters, data):

    """Apply useful filters on the data
    """
    opers = {'<' : operator.lt,
             '=' : operator.eq,
             '>' : operator.gt}
    
    for filter in filters:
        op= filter.split()[1]
        key = filter.split()[0]
        val = float(filter.split()[2])
        data = data[opers[op](data[key], val)]
    
    return data

if __name__ == '__main__':
    
    if len(sys.argv) < 4:
        print('Usage: datapath parse/load fuel_cool')
        sys.exit()

    args = sys.argv[1:]
    # get datapaths
    path = args[0]
    trainpath = path + args[2] + '_train'
    testpath = path + args[2] + '_test'
    
    if args[1] == 'parse':
        save_store_data(trainpath + '.zip', trainpath + '.csv')    
        save_store_data(testpath + '.zip', testpath + '.csv')
    
    # load training and testing sets
    traindata = load_from_csv('{0}.csv'.format(trainpath))
    testdata = load_from_csv('{0}.csv'.format(testpath))
    # test the interpolator
    func = interpolate_grid(traindata)
    test_interpolator(testdata, func, args[2])
