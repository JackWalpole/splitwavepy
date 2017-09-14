"""
Some plotting routines
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt

def plot_data(data):
    """
    Plot trace data stored in (2,n) numpy array
    """
    from matplotlib import gridspec
    fig = plt.figure(figsize=(12, 3)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    # trace data
    ax0 = plt.subplot(gs[0])
    ax0.plot(data.T)
    # particle  motion
    lim = abs(data).max() * 1.1
    # the polar axis:
    # ax_polar = plt.subplot(gs[1], polar=True, frameon=False)
    # ax_polar.set_rmax(lim)
    # ax_polar.patch.set_facecolor(111)
    # ax_polar.get_xaxis.set_visible(False)
    # ax_polar.grid(True)
    # the data
    ax1 = plt.subplot(gs[1])
    # ax1.patch.set_alpha(0)
    ax1.plot(data[1],data[0])
    ax1.set_xlim([-lim,lim])
    ax1.set_ylim([-lim,lim])
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    # show
    plt.show()
    
def plot_surf(X,Y,Z,cmap='viridis'):
    plt.contourf(X,Y,Z,cmap=cmap)
    plt.show()