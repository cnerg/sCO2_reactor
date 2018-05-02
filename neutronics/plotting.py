import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy as np
import pandas as pd

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns


axis_labels = {'core_r' : 'Core R [cm]',
               'cool_r' : 'Coolant Channel R [cm]',
               'power' : 'Power [kW]',
               'PD' : 'Pitch/Cool. D. [-]',
               'enrich' : 'enrich [-]',
               'keff' : 'EOL keff [-]',
               'ave_E' : 'Average Energy EOL [MeV]',
               'mass' : 'Fuel mass [kg]'
              }

data = pd.read_csv("depl_results.csv")
#data = filter_data([('keff', 'great', 1.0)], data)

def plot_results(ind, dep, colorplot=[]):
    """Generate Plots
    """
    # plot
    fig = plt.figure()
    if len(colorplot) > 0:
        plt.scatter(ind, dep, c=colorplot, s=6,
                cmap=plt.cm.get_cmap('plasma', len(set(colorplot))))
        plt.colorbar()
    else:
        plt.scatter(ind, dep, s=6)
    
    plt.show()

    return plt

def property_plot_matrix(data, params, colorplot=None):
        """
        """
        
        dim = len(params)

        plots = []
        for i in range(dim):
            row = []
            for j in range(dim):
                row.append((params[i], params[j]))
            plots.append(row)

        pass_idx = []
        # create list of "skippable" repeat plots
        for rowdx in range(dim):
            basedx = rowdx * dim 
            for i in range(rowdx):
                pass_idx.append(basedx + i)


        fig = plt.figure(figsize=(8, 6))
        for xidx, row in enumerate(plots):
            for yidx, plot in enumerate(row): 
                pass_id = dim*yidx + xidx
                if pass_id not in pass_idx:
                    ax = fig.add_subplot(dim, dim, pass_id+1)
                    xkey = plot[0] 
                    ykey = plot[1]
                    x = data[xkey]
                    y = data[ykey]
                    ax.tick_params(which="both", size=0.25)
                    if colorplot:
                        ax.scatter(x, y, s=(3/dim), c=data[colorplot], cmap=plt.cm.get_cmap('plasma', len(set(data['keff']))))
                    else:
                        ax.scatter(x, y, s=(3/dim), cmap=plt.cm.get_cmap('plasma', len(set(data['keff']))))
                    if ykey == 'power':
                        ax.ticklabel_format(axis="y", style="sci", fontsize=8, scilimits=(0, 0))
                    if xkey == 'power':
                        ax.ticklabel_format(axis="x", style="sci", fontsize=8, scilimits=(0, 0))
                    if yidx == 0:
                        plt.title(axis_labels[xkey], fontsize=8)
                    if xidx == yidx:
                        plt.ylabel(axis_labels[ykey], fontsize=8)

                    if xidx != yidx:
                        ax.xaxis.set_visible(False)
                        ax.yaxis.set_visible(False)
        fig.savefig('sampling_plot_matrix.png', dpi=500)
        
        
        return plt

def property_plot_row(data, params, dep, colorplot=None):
        """
        """
        
        dim = len(params)
        basesize = 6
        fig = plt.figure(figsize=(20, 10))
        for xidx, plot in enumerate(params):
            ax = fig.add_subplot(1, dim, xidx+1)
            xkey = plot
            ykey = dep
            x = data[xkey]
            y = data[ykey]
            ax.tick_params(which="both")
            
            if colorplot:
                im = ax.scatter(x, y, s=(basesize/dim), c=data[colorplot],
                            cmap=plt.cm.get_cmap('plasma', 
                            len(set(data[colorplot]))))
            else:
                ax.scatter(x, y, s=(basesize/dim))
            
            ax.set_xlabel(axis_labels[xkey])
            ax.set_ylabel(axis_labels[ykey])    
            plt.title(xkey, fontsize=12)
            
            if xidx != 0:
                ax.yaxis.set_visible(False)
            
        
        fig.savefig('sampling_plot_row.png', dpi=500)
        
        if colorplot:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.007, 0.7])
            fig.colorbar(im, cax=cbar_ax, label=axis_labels[colorplot])

        plt.show() 
        
        return plt

def surf_plot(data, ind1, ind2, dep, colorplot=None):
    """
    """

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X = data[ind1]
    Y = data[ind2]
    Z = data[dep]
    scaled_size = np.log(np.multiply(data['keff'], 10))*10

    if colorplot:
        # Plot the surface.
        p = ax.scatter(X,Y,Z, s=scaled_size, c=data[colorplot],
                cmap=plt.cm.get_cmap('viridis',
            len(data[colorplot])))
        plt.colorbar(p, ax=ax, label=axis_labels[colorplot])
    else:
        ax.scatter(X, Y, Z, c=Z)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    labelsize = 5
    ax.set_xlabel(axis_labels[ind1])
    ax.set_ylabel(axis_labels[ind2])
    ax.set_zlabel(axis_labels[dep])
    
    ax.xaxis.label.set_size(labelsize)
    ax.yaxis.label.set_size(labelsize)
    ax.zaxis.label.set_size(labelsize)
    
    plt.show()
    plt.savefig('surface_plot.png', dpi=500)


#plt = property_plot_matrix(data, ('core_r', 'PD', 'enrich', 'mass', 'keff'))
plt = property_plot_row(data, ['core_r', 'PD', 'enrich', 'cool_r'], 'keff', 'mass')
#surf_plot(data, 'core_r', 'enrich', 'PD', 'keff')
