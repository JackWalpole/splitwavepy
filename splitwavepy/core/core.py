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

import numpy as np
import math
from scipy import signal


def lag(data,nsamps):
    """
    Lag t1 nsamps/2 to the left and
    lag t2 nsamps/2 to the right.
    nsamps must be even.
    If nsamps is negative t1 shifted to the right and t2 to the left.
    This process truncates trace length by nsamps and preserves centrality.
    Therefore windowing must be used after this process 
    to ensure even trace lengths when measuring splitting.
    """
    if nsamps == 0:
        return data
    elif nsamps%2 != 0:
        raise Exception('nsamps must be even')
    
    t1 = data[0,:]
    t2 = data[1,:]
    
    if nsamps > 0:
        # positive shift
        return np.vstack((t1[nsamps:], t2[:-nsamps]))
    else:
        # negative shift
        return np.vstack((t1[:nsamps], t2[-nsamps:]))

        
def rotate(data,degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates from x to y axis
       e.g. N to E if row 0 is N cmp and row1 is E cmp"""
    ang = np.deg2rad(degrees)
    rot = np.array([[np.cos(ang),-np.sin(ang)],
                    [np.sin(ang), np.cos(ang)]])
    return np.dot(rot,data)

# def rotate_and_lag(data,degrees,nsamps):
#     return lag(rotate(data,degrees), nsamps)

def split(data,degrees,nsamps):
    """Apply forward splitting and rotate back"""
    data = rotate(data,-degrees)
    data = lag(data,nsamps)
    data = rotate(data,degrees)
    return data

def unsplit(data,degrees,nsamps):
    """Apply inverse splitting and rotate back"""
    return split(data,degrees,-nsamps)
    
def chop(data,nsamps,centre=None,tukey=None):
    """Chop trace, or traces, about *centre* sample (defaults to middle sample)
       width must be odd to maintain definite centre"""
    
    if data.ndim == 1:
        length = data.shape[0]
    else:
        length = data.shape[1]
    
    if nsamps%2 != 1:
        raise Exception('width must be odd')
    
    if centre == None:
        # by defualt put centre at middle sample
        centre = int(length/2)
    else:
        centre = int(centre)
        
    if nsamps > length:
        raise Exception('chop width is greater than trace length')
    
    hw = int(nsamps/2)    
    t0 = centre - hw
    t1 = centre + hw
    
    if t0 < 0:
        raise Exception('chop starts before trace data')
    elif t1 > length:
        raise Exception('chop ends after trace data')
        
    if tukey is not None:
        tukey = signal.tukey(nsamps,alpha=tukey)
    else:
        tukey = 1.
        
    if data.ndim == 1:
        return data[t0:t1+1] * tukey
    else: 
        return data[:,t0:t1+1] * tukey
    
def pca(data):
    """
    Principal component analysis
    Returns direction of strongest component in degrees anti-clockwise from x
    """
    w,v = np.linalg.eig(np.cov(data))
    m = np.argmax(w)
    return np.rad2deg(np.arctan2(v[1,m],v[0,m]))
    
def snr(data):
    """
    Returns signal to noise ratio assuming signal on trace1 and noise on trace2
    Calculates on each trace by sum of squares and then takes the ratio
    """
    signal = np.sum(data[0,:]**2)
    noise = np.sum(data[1,:]**2)
    return signal / noise
    
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

def synth(pol=0,fast=0,lag=0,noise=0.05,nsamps=501,width=16.0):
    """return ricker wavelet synthetic data"""
    ricker = signal.ricker(int(nsamps), width)
    data = np.vstack((ricker,np.zeros(ricker.shape)))
    # white noise convolved with a gaussian wavelet
    noise = np.random.normal(0,noise,data.shape)
    std = width/4
    norm = 1/(std*np.sqrt(2*np.pi))
    gauss = norm * signal.gaussian(int(nsamps),std)
    noise[0] = np.convolve(noise[0],gauss,'same')
    noise[1] = np.convolve(noise[1],gauss,'same')
    data = data + noise
    data = rotate(data,pol)
    return split(data,fast,lag)
    
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
