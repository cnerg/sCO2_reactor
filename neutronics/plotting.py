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
# apply linear or log shift to data
ops = {'lin' : lambda x: x,
       'log' : lambda x: np.log(x)
      }

data = pd.read_csv("depl_results.csv")
#data = filter_data([('keff', 'great', 1.0)], data)

def plot_results(ind, dep, labels):
    """Generate Plots
    """
    fig = plt.figure()
     
    for x, y, label in zip(ind, dep, labels):
        # plot
        plt.plot(x, y, label=label)
    
    plt.legend(labels) 
    plt.show()

    return plt

def property_plot_matrix(data, params, colorplot=None):
        """Produce a matrix property plot, plotting every data category against
        itself and the other categories

        Arguments:
        ---------- 
            data (array):  panda array of the all data
            params (tuple): tuple of keys for desired plots
            colorplot (str): optional str to include colorbar
        Returns:
        --------
            plt (matplotlib inst.): plot object
        """
        
        dim = len(params)

        plots = []
        for i in range(dim):
            row = []
            for j in range(dim):
                row.append((params[i], params[j]))
            plots.append(row)

        pass_idx = []
        # create list of "skippable" repeat plots, forming RU 'matrix' of plots
        for rowdx in range(dim):
            basedx = rowdx * dim 
            for i in range(rowdx):
                pass_idx.append(basedx + i)


        fig = plt.figure(figsize=(8, 6))
        # loop throug each plot, creating it and assigning it to the grid
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
                        ax.scatter(x, y, s=(3/dim), 
                                   c=data[colorplot], 
                                   cmap=plt.cm.get_cmap('plasma', 
                                     len(set(data[colorplot]))))
                    else:
                        ax.scatter(x, y, s=(3/dim))
                    if ykey == 'power':
                        ax.ticklabel_format(axis="y", style="sci", 
                                            fontsize=8, scilimits=(0, 0))
                    if xkey == 'power':
                        ax.ticklabel_format(axis="x", style="sci", 
                                            fontsize=8, scilimits=(0, 0))
                    if yidx == 0:
                        plt.title(axis_labels[xkey], fontsize=8)
                    if xidx == yidx:
                        plt.ylabel(axis_labels[ykey], fontsize=8)
                    # only display axis labels on center diagonal, saves space
                    # on matrix
                    if xidx != yidx:
                        ax.xaxis.set_visible(False)
                        ax.yaxis.set_visible(False)
        fig.savefig('sampling_plot_matrix.png', dpi=500)
        
        
        return plt

def property_plot_row(data, params, dep, loglin, colorplot=None):
        """
        """
        
        dim = len(params)
        basesize = 6
        fig = plt.figure(figsize=(12, 8))
        for xidx, plot in enumerate(params):
            ax = fig.add_subplot(1, dim, xidx+1)
            xkey = plot
            ykey = dep
            x = ops[loglin[0]](data[xkey])
            y = ops[loglin[1]](data[ykey])
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

if __name__=='__main__':
    #plt = property_plot_matrix(data, ('core_r', 'PD', 'enrich', 'mass', 'keff'))
    plt = property_plot_row(data, ['core_r', 'PD'], 'mass', ('lin', 'lin'),
    'keff')
    #surf_plot(data, 'core_r', 'enrich', 'PD', 'keff')
