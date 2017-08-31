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

from . import plotting as p
from . import eigval

# Classes

class Pair:
    """
    The Pair is a class to store two traces in the x and y directions.
    Methods are included to facilitate analysis on this Pair of traces.
    If data is not provided on initiation will return a ricker wavelet with noise.
    Usage: Pair()     => create Pair of synthetic data
           Pair(data) => creates Pair from two traces stored as rows in numpy array data
           Pair(x,y) => creates Pair from two traces stored in numpy arrays x and y.
    """
    def __init__(self,*args):
        
        if len(args) == 0:                      
            self.data = synth()            
        elif len(args) == 1:            
            self.data = args[0]       
        elif len(args) == 2:            
            self.data = np.vstack((args[0],args[1]))     
        else: 
            raise Exception('Unexpected number of arguments')
                    
        # some sanity checks
        if self.data.ndim != 2:
            raise Exception('data must be two dimensional')
        if self.data.shape[0] != 2:
            raise Exception('data must contain two traces in two rows')
        if self.data.shape[1]%2 == 0:
            raise Exception('traces must have odd number of samples')

        
    # methods
    def plot(self):
        p.plot_data(self.data)
    
    def split(self,degrees,nsamps):
        self.data = split(self.data,degrees,nsamps)
        return self
    
    def unsplit(self,degrees,nsamps):
        self.data = unsplit(self.data,degrees,nsamps)
        return self
        
    def rotate(self,degrees):
        self.data = rotate(self.data,degrees)
        return self
        
    def lag(self,nsamps):
        self.data = lag(self.data,nsamps)
        return self
        
    def window(self,width):
        self.data = window(self.data,width)
        return self
        
    def grideigval(self, maxshift=None, window=None, stepang=None, stepshift=None):
        return eigval.grideigval(self.data)



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
       rotates from x to y axis"""
    ang = np.deg2rad(degrees)
    rot = np.array([[np.cos(ang),-np.sin(ang)],
                    [np.sin(ang), np.cos(ang)]])
    return np.dot(rot,data)

def rotate_and_lag(data,degrees,nsamps):
    return lag(rotate(data,degrees), nsamps)

def split(data,degrees,nsamps):
    """Apply forward splitting and rotate back"""
    return rotate( rotate_and_lag(data,degrees,nsamps), -degrees)

def unsplit(data,degrees,nsamps):
    """Apply inverse splitting and rotate back"""
    return rotate( rotate_and_lag(data,degrees,-nsamps), -degrees)
    
def window(data,width):
    """Window trace, or traces, about centre sample
       width must be odd to maintain definite centre"""
    
    if data.ndim == 1:
        length = data.shape[0]
    else:
        length = data.shape[1]
    
    if width%2 != 1:
        raise Exception('width must be odd')
    
    centre = math.floor(length/2)
        
    if width > length:
        raise Exception('window width is greater than trace length')
        
    t0 = centre - math.floor(width/2)
    t1 = centre + math.ceil(width/2)
    
    if t0 < 0:
        raise Exception('window starts before trace data')
    elif t1 > length:
        raise Exception('window ends after trace data')
     
    if data.ndim == 1:
        return data[t0:t1]
    else: 
        return data[:,t0:t1]
    
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

def synth(srcpol=0,fast=0,lag=0,noise=0.005,nsamps=501,width=16.0):
    """return ricker wavelet synthetic data"""
    ricker = signal.ricker(int(nsamps), width)
    data = np.vstack((ricker,np.zeros(ricker.shape)))
    data = data + np.random.normal(0,noise,data.shape)
    data = rotate(data,srcpol)
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
