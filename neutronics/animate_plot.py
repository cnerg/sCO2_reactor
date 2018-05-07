import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
        p = ax.scatter(X,Y,Z, s=np.log(data['mass'])*50, c=data[colorplot], cmap=plt.cm.get_cmap('viridis',
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
    return fig,

# Animate
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=True, repeat=True)
# Save
anim.save('basic_animation.mp4', fps=10, bitrate=-1,
        dpi=500,extra_args=['-vcodec', 'libx264'], writer='ffmpeg')
