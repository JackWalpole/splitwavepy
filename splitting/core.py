"""
Low level routines for dealing with traces in numpy arrays.
Works on arrays sample by sample and need not know 
anything about the time of a sample interval.
Assumes data of interest is at the centre of the array.
Trace must have odd number of samples so there is a definite centre.
"""

from __future__ import print_function

import numpy as np
import math
from scipy import signal

from . import plotting as p

# Classes

class Pair:
    """
    The Pair is the two traces on which analysis is performed
    it doesn't need to be used but supports methods which may
    make the syntax easier to use
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
        self.plot = p.plot_pair(self.data)
    
    def split(self,degrees,nsamps):
        self.data = split(self.data,degrees,nsamps)
    
    def unsplit(self,degrees,nsamps):
        self.data = unsplit(self.data,degrees,nsamps)
        
    def rotate(self,degrees):
        self.data = rotate(self.data,degrees)
        
    def lag(self,nsamps):
        self.data = lag(self.data,nsamps)
        

def synth(noise=True,fast=0,lag=0):
    """return ricker wavelet synthetic data"""
    ricker = signal.ricker(501, 16.0)
    pair = np.vstack((ricker,np.zeros(ricker.shape)))
    if noise is True:
        noise = np.random.normal(0,.005,pair.shape)
        pair = pair + noise
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
    
# def hann_window(pair,width):
#     """Apply Hanning taper
#        of width number of samples, truncates traces"""
#     data = window(pair,width)
#     return np.hanning(width) * data

# def eigcov(pair):
#     return np.sort(np.linalg.eigvals(np.cov(pair,rowvar=True)))
# #     return np.sort(np.linalg.eigvals(np.cov(pair)))
#
# def grideigcov(pair,maxshift,window=None, stepang=None,stepshift=None):
#
#     # set some defaults
#     if stepshift is None:
#         stepshift = 2 * int(np.max([1,maxshift/40]))
#     if stepang is None:
#         stepang = 2
#     if window is None:
#         # by default whatevers smaller,
#         # half trace length or 10 * max shift
#         window = int(np.min([pair.shape[1] * 0.5,maxshift * 10]))
#
#     deg, lag = np.meshgrid(np.arange(0,180,stepang),
#                              np.arange(0,maxshift,stepshift).astype(int))
#
#     shape = deg.shape
#     lam1 = np.zeros(shape)
#     lam2 = np.zeros(shape)
#     for ii in np.arange(shape[1]):
#         temp = rotate(pair,deg[0,ii]+90)
#         for jj in np.arange(shape[0]):
#             temp2 = shift(temp,lag[jj,ii])
#             temp3 = taper(temp2,window)
#             lam2[jj,ii], lam1[jj,ii] = eigcov(temp3)
#     return deg, lag, lam1, lam2
#
# def xcorr(pair):
#     return np.correlate(pair[0,:],pair[1,:])[0]
#
# def gridxcorr(pair,maxshift,window=None, stepang=None,stepshift=None):
#
#     # set some defaults
#     if stepshift is None:
#         stepshift = int(np.max([1,maxshift/40]))
#     if stepang is None:
#         stepang = 2
#     if window is None:
#         # by default whatevers smaller,
#         # half trace length or 10 * max shift
#         window = int(np.min([pair.shape[1] * 0.5,maxshift * 10]))
#
#     deg, shift = np.meshgrid(np.arange(0,180,stepang),
#                              np.arange(0,maxshift,stepshift).astype(int))
#
#     shape = deg.shape
#     xc = np.zeros(shape)
#     for ii in np.arange(shape[1]):
#         temp = rotate(pair,deg[0,ii]+90)
#         for jj in np.arange(shape[0]):
#             temp2 = shift(temp,shift[jj,ii])
#             temp3 = taper(temp2,window)
#             xc[jj,ii] = xcorr(temp3)
#     return deg, shift, xc
    