"""
Low level routines for dealing with traces in numpy arrays.
Works on arrays sample by sample and need not know 
anything about the time of a sample interval.
Assumes data of interest is at the centre of the array.
Trace must have odd number of samples so there is a definite centre.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .window import Window

import numpy as np
from scipy import signal

def lag(x,y,nsamps):
    """
    Lag x nsamps/2 to the left and
    lag y nsamps/2 to the right.
    nsamps must be even.
    If nsamps is negative x shifted to the right and y to the left.
    This process truncates trace length by nsamps and preserves centrality.
    Therefore windowing must be used after this process 
    to ensure even trace lengths when measuring splitting.
    """
    if nsamps == 0:
        return x,y
    elif nsamps%2 != 0:
        raise Exception('nsamps must be even')

    if nsamps > 0:
        # positive shift
        return x[nsamps:], y[:-nsamps]
    else:
        # negative shift
        return x[:nsamps], y[-nsamps:]

        
def rotate(x,y,degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates from x to y axis
       e.g. N to E if row 0 is N cmp and row1 is E cmp"""
    ang = np.deg2rad(degrees)
    rot = np.array([[np.cos(ang),np.sin(ang)],
                    [-np.sin(ang), np.cos(ang)]])
    xy = np.dot(rot,np.vstack((x,y)))
    return xy[0], xy[1]

# def rotate_and_lag(data,degrees,nsamps):
#     return lag(rotate(data,degrees), nsamps)

def split(x,y,degrees,nsamps):
    """Apply forward splitting and rotate back"""
    if nsamps == 0:
        return x,y
    x,y = rotate(x,y,degrees)
    x,y = lag(x,y,nsamps)
    x,y = rotate(x,y,-degrees)
    return x,y

def unsplit(x,y,degrees,nsamps):
    """Apply inverse splitting and rotate back"""
    return split(x,y,degrees,-nsamps)

def chop(*args,**kwargs):
    """Chop trace, or traces, using window"""
    
    if ('window' in kwargs):
        window = kwargs['window']
    
    if not isinstance(window,Window):
        raise Exception('window must be a Window')
    
    length = args[0].size
          
    if window.width > length:
        raise Exception('window width is greater than trace length')
    
    centre = int(length/2) + window.offset
    hw = int(window.width/2)    
    t0 = centre - hw
    t1 = centre + hw
    
    if t0 < 0:
        raise Exception('chop starts before trace data')
    elif t1 > length:
        raise Exception('chop ends after trace data')
        
    if window.tukey is not None:
        tukey = signal.tukey(window.width,alpha=window.tukey)
    else:
        tukey = 1.
    
    if len(args)==1:    
        return args[0][t0:t1+1] * tukey
    elif len(args)==2:
        return args[0][t0:t1+1] * tukey, args[1][t0:t1+1] * tukey
    elif len(args)==3:
        return args[0][t0:t1+1] * tukey, args[1][t0:t1+1] * tukey, args[2][t0:t1+1] * tukey

    
def pca(x,y):
    """
    Principal component analysis
    Returns direction of strongest component in degrees anti-clockwise from x
    """
    w,v = np.linalg.eig(np.cov(np.vstack((x,y))))
    m = np.argmax(w)
    return np.rad2deg(np.arctan2(v[1,m],v[0,m]))
    
# def snr(data):
#     """
#     Returns signal to noise ratio assuming signal on trace1 and noise on trace2
#     Calculates on each trace by sum of squares and then takes the ratio
#     """
#     signal = np.sum(data[0,:]**2)
#     noise = np.sum(data[1,:]**2)
#     return signal / noise
    
def snrRH(data):
    """
    Returns signal to noise ratio assuming signal on trace1 and noise on trace2
    Uses the method of Restivo and Helffrich (1999):
    peak amplitude on trace1 / 2*std trace2
    """
    signal = np.max(data[0,:])
    noise = 2 * np.std(data[1,:])
    return signal / noise

# Useful bits and pieces
    
def noise(size,amp,smooth):
    """gaussian noise convolved with a gaussian wavelet"""
    norm = 1/(smooth*np.sqrt(2*np.pi))
    gauss = norm * signal.gaussian(size,smooth)
    n = np.random.normal(0,amp,size)
    return np.convolve(n,gauss,'same')  
    
def min_idx(vals):
    """
    return indices of min value in vals grid
    """
    return np.unravel_index(np.argmin(vals),vals.shape)
    
def max_idx(vals):
    """
    return indice of max value in vals grid
    """
    return np.unravel_index(np.argmax(vals),vals.shape)
