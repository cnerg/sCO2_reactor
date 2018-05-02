import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import animation

import numpy as np
import pandas as pd

from parse_outputs import filter_data, plot_results
import neutronic_sweeps as ns


data = pd.read_csv("depl_results.csv")
data = filter_data([('keff', 'great', 1.0)], data)

fig = plt.figure()
ax = fig.gca(projection='3d')

axrotate = plt.axes([0.15, 0.02, 0.5, 0.03], facecolor='lightgoldenrodyellow')


def surf_plot(data, ind1, ind2, dep, colorplot=None):
    """
    """
    axis_labels = {'core_r' : 'Core r [cm]',
                   'power' : 'Power [kW]',
                   'PD' : 'Pitch/Cool. D. [-]',
                   'enrich' : 'enrich [-]',
                   'keff' : 'EOL keff [-]',
                   'ave_E' : 'Average Energy EOL [MeV]',
                   'mass' : 'Fuel mass [kg]'
                   }


    X = data[ind1]
    Y = data[ind2]
    Z = data[dep]
    if colorplot:
        # Plot the surface.
        p = ax.scatter(X,Y,Z, s=np.log(data['keff'])*500, c=data[colorplot], cmap=plt.cm.get_cmap('viridis',
            len(data[colorplot])))
        plt.colorbar(p, ax=ax, label=axis_labels[colorplot])
    else:
        ax.scatter(X, Y, Z, c=Z)

    # Customize the z axis.
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    ax.set_xlabel(axis_labels[ind1])
    ax.set_ylabel(axis_labels[ind2])
    ax.set_zlabel(axis_labels[dep])


def init():
    surf_plot(data, 'core_r', 'enrich', 'PD', 'keff')
    return fig,

def animate(i):
    ax.view_init(elev=10., azim=i)
    fig.canvas.draw_idle()
    return fig,

init()
srotate = Slider(axrotate, 'Rotate', 0, 720, valinit=0)
srotate.valtext.set_visible(False)
srotate.on_changed(animate)
plt.show()
plt.savefig('slider.png', dpi=700)
