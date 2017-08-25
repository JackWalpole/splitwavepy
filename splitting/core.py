"""
Low level routines for dealing with traces in numpy arrays.
Works on arrays sample by sample and need not know 
anything about the time of a sample interval.
Assumes data of interest is at the centre of the array.
Trace must have odd number of samples so there is a definite centre.
"""

import numpy as np
import math
from scipy import signal

from . import plotting as p
from . import eigval

# Classes

class Pair:
    """
    The Pair is a convenience class in which to store two traces x and y.
    Methods are included to facilitate analysis on this Pair of traces.
    If x and y data not provided will generate a ricker wavelet with noise.
    """
    def __init__(self,x=None,y=None):
        """x and y are traces in the form of numpy arrays"""
        
        if x is None or y is None:
            
            self.data = synth()
        
        else:
            
            # some checks
            if x.ndim != 1:
                raise Exception('traces must be one dimensional')
            if x.shape != y.shape:
                raise Exception('x and y must be the same shape')
            if x.shape[0]%2 == 0:
                raise Exception('traces must have odd number of samples')
            # combine data into pair
            self.data = np.vstack((x,y))
        
    # methods
    def plot(self):
        p.plot_pair(self.data)
    
    def split(self,degrees,nsamps):
        self.data = split(self.data,degrees,nsamps)
    
    def unsplit(self,degrees,nsamps):
        self.data = unsplit(self.data,degrees,nsamps)
        
    def rotate(self,degrees):
        self.data = rotate(self.data,degrees)
        
    def lag(self,nsamps):
        self.data = lag(self.data,nsamps)
        
    def window(self,width):
        self.data = window(slef.data,width)
        
    def grideigval(self, maxshift=None, window=None, stepang=None, stepshift=None):
        return eigval.grideigvalcov(self.data)
        

def synth(fast=0,lag=0,noise=0.005,nsamps=501,width=16.0):
    """return ricker wavelet synthetic data"""
    ricker = signal.ricker(int(nsamps), width)
    pair = np.vstack((ricker,np.zeros(ricker.shape)))
    pair = pair + np.random.normal(0,noise,pair.shape)
    return split(pair,fast,lag)

def lag(pair,nsamps):
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
        return pair
    elif nsamps%2 != 0:
        raise Exception('nsamps must be even')
    
    t1 = pair[0,:]
    t2 = pair[1,:]
    
    if nsamps > 0:
        # positive shift
        return np.vstack((t1[nsamps:], t2[:-nsamps]))
    else:
        # negative shift
        return np.vstack((t1[:nsamps], t2[-nsamps:]))

        
def rotate(pair,degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates clockwise"""
    ang = np.deg2rad(degrees)
    rot = np.array([[np.cos(ang),-np.sin(ang)],
                    [np.sin(ang), np.cos(ang)]])
    return np.dot(rot,pair)

def rotate_and_lag(pair,degrees,nsamps):
    return lag(rotate(pair,degrees), nsamps)

def split(pair,degrees,nsamps):
    """Apply forward splitting and rotate back"""
    return rotate( rotate_and_lag(pair,degrees,nsamps), -degrees)

def unsplit(pair,degrees,nsamps):
    """Apply inverse splitting and rotate back"""
    return rotate( rotate_and_lag(pair,degrees,-nsamps), -degrees)
    
def window(pair,width):
    """Window traces about centre sample
       width must be odd to maintain definite centre"""
    
    if width%2 != 1:
        raise Exception('width must be odd')
    
    centre = math.floor(pair.shape[1]/2)
        
    if width > pair.shape[1]:
        raise Exception('window width is greater than trace length')
        
    t0 = centre - math.floor(width/2)
    t1 = centre + math.ceil(width/2)
    
    if t0 < 0:
        raise Exception('window starts before trace data')
    elif t1 > pair.shape[1]:
        raise Exception('window ends after trace data')
        
    return pair[:,t0:t1]

# Useful bits and pieces
    
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
