"""
Routines to work with 3-component data
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .window import Window
from . import geom
from . import core

import numpy as np
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
    if nsamps%2 != 0:
        raise Exception('nsamps must be even')
    
    if nsamps == 0:
        return data
    else:
        return np.vstack((core.lag(data[0:2],nsamps),data[2][int(nsamps/2):-int(nsamps/2)]))
        
def rotate(data,degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates from x to y axis
       e.g. N to E if row 0 is N cmp and row1 is E cmp"""
    return np.vstack((core.rotate(data[0:2],degrees),data[2]))

def split(data,degrees,nsamps):
    """Apply forward splitting and rotate back"""
    data = rotate(data,-degrees)
    data = lag(data,nsamps)
    data = rotate(data,degrees)
    return data

def unsplit(data,degrees,nsamps):
    """Apply inverse splitting and rotate back"""
    return split(data,degrees,-nsamps)
    
def chop(data,window):
    """Chop trace, or traces, using window"""
    return core.chop(data,window)

# def pca(data):
#     """
#     Principal component analysis
#     Returns direction of strongest component in degrees anti-clockwise from x
#     """
#     w,v = np.linalg.eig(np.cov(data))
#     m = np.argmax(w)
#     return np.rad2deg(np.arctan2(v[1,m],v[0,m]))

def synth(pol=0,fast=0,lag=0,noise=0.05,nsamps=501,width=16.0,**kwargs):
    """return ricker wavelet synthetic data"""
    ricker = signal.ricker(int(nsamps), width)
    data = np.vstack((ricker,np.zeros((2,ricker.shape[0]))))
    # gaussian noise convolved with a gaussian wavelet
    noise = np.random.normal(0,noise,data.shape)
    std = width/4
    norm = 1/(std*np.sqrt(2*np.pi))
    gauss = norm * signal.gaussian(int(nsamps),std)
    noise[0] = np.convolve(noise[0],gauss,'same')
    noise[1] = np.convolve(noise[1],gauss,'same')
    noise[2] = np.convolve(noise[2],gauss,'same')
    data = data + noise
    data = rotate(data,pol)
    data = split(data,fast,lag)
    return data
    