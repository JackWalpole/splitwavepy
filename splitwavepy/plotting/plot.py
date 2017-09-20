"""
Some plotting routines
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# matplotlib stuff
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D

import numpy as np


def trace(*args,**kwargs):
    """Return axis with trace data.
    
    kwargs:
    - time
    - ax
    """
    
    if 'time' not in kwargs:
        kwargs['time'] = np.arange(args[0].size)   
    
    if 'ax' not in kwargs:        
        # initiate axis  
        kwargs['ax'] = plt.subplot(111)    
    ax = kwargs['ax']
        
    # plot data
    for ii in range(len(args)):
        ax.plot(kwargs['time'],args[ii])
        
    # set limit
    lim = np.abs(args).max() * 1.1
    ax.set_ylim([-lim,lim])    
        
    # set label
    if 'units' in kwargs:
        units = kwargs['units']
    else:
        units = 's'    
    ax.set_xlabel('Time (' + units +')')  
    
    # plot window markers
    if 'window' in kwargs:
        nsamps = args[0].size
        wbeg = kwargs['window'].start(nsamps)*kwargs['time'][1]
        wend = kwargs['window'].end(nsamps)*kwargs['time'][1]
        ax.axvline(wbeg,linewidth=1,color='k')
        ax.axvline(wend,linewidth=1,color='k')        
    
    return ax
    

def particle(*args,**kwargs):
    """Return axis with particle motion data
    
    kwargs:
    - ax
    """
    if not ('labels' in kwargs):
        kwargs['labels'] = ['x','y','z']
    
    # 2D particle motion
    if len(args) == 2:
        if 'ax' not in kwargs:         
            kwargs['ax'] = plt.subplot(111)
        ax = kwargs['ax']
        ax.plot(args[1],args[0])
        lim = np.abs(args).max() * 1.1
        ax.set_aspect('equal')
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
        ax.set_xlabel(kwargs['labels'][1])
        ax.set_ylabel(kwargs['labels'][0])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.grid()
        return ax
    
    # 3D particle motion
    if len(args) == 3:
        if 'ax' not in kwargs:
            kwargs['ax'] = plt.subplot(111,projection='3d')
        ax = kwargs['ax']
        ax.plot(args[0],args[1],args[2])
        lim = np.abs(args).max() * 1.1
        ax.set_aspect('equal')
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
        ax.set_zlim([-lim,lim])
        ax.plot(args[0],args[1],-lim,zdir='z',alpha=0.3,color='g')
        ax.plot(args[0],args[2],lim,zdir='y',alpha=0.3,color='g')
        ax.plot(args[1],args[2],-lim,zdir='x',alpha=0.3,color='g')
        ax.set_xlabel(kwargs['labels'][0])
        ax.set_ylabel(kwargs['labels'][1])
        ax.set_zlabel(kwargs['labels'][2])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.axes.zaxis.set_ticklabels([])
        return ax
        
        
        

# def surface(*args,**kwargs):
#     """Plot 2-D surface data"""
#
#
# def plot_data(data):
#     """
#     Plot trace data stored in (2,n) numpy array
#     """
#     from matplotlib import gridspec
#     fig = plt.figure(figsize=(12, 3))
#     gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
#     # trace data
#     ax0 = plt.subplot(gs[0])
#     ax0.plot(data.T)
#     # particle  motion
#     lim = abs(data).max() * 1.1
#     # the polar axis:
#     # ax_polar = plt.subplot(gs[1], polar=True, frameon=False)
#     # ax_polar.set_rmax(lim)
#     # ax_polar.patch.set_facecolor(111)
#     # ax_polar.get_xaxis.set_visible(False)
#     # ax_polar.grid(True)
#     # the data
#     ax1 = plt.subplot(gs[1])
#     # ax1.patch.set_alpha(0)
#     ax1.plot(data[1],data[0])
#     ax1.set_xlim([-lim,lim])
#     ax1.set_ylim([-lim,lim])
#     ax1.axes.get_xaxis().set_visible(False)
#     ax1.axes.get_yaxis().set_visible(False)
#     # show
#     plt.show()
#
# def plot_surf(X,Y,Z,cmap='viridis'):
#     plt.contourf(X,Y,Z,cmap=cmap)
#     plt.show()



