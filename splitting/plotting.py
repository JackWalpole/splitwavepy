"""
Some plotting routines
"""

import matplotlib.pyplot as plt

def plot_data(data):
    """
    Plot trace data stored in (2,n) numpy array
    """
    from matplotlib import gridspec
    fig = plt.figure(figsize=(12, 3)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    ax0 = plt.subplot(gs[0])
    ax0.plot(data.T)
    ax1 = plt.subplot(gs[1])
    ax1.plot(data[0,:],data[1,:])
    plt.axis('equal')
    plt.show()
    
def plot_surf(X,Y,Z,cmap='viridis'):
    plt.contourf(X,Y,Z,cmap=cmap)
    plt.show()