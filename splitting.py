#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal

def shift(pair,shift):
    """shift t1 shift/2 samples to the left and
       shift t2 shift/2 samples to the right,
       shift must be even number of samples
       this process truncates trace length"""
    if shift == 0:
        return pair
    elif shift%2 != 0:
        raise Exception('shift must be even number of samples')
    
    if shift > 0:
        t1 = pair[0,:]
        t2 = pair[1,:]
    else:
        t1 = pair[1,0]
        t2 = pair[0,:]

    return np.vstack((t1[math.ceil(shift/2):-math.floor(shift/2)], t2[:-shift]))
        


def rotate(pair,degrees):
    """t1 is x-axis and t2 is y-axis,
       rotates clockwise"""
    ang = np.deg2rad(degrees)
    rot = np.array([[np.cos(ang),-np.sin(ang)],
                    [np.sin(ang), np.cos(ang)]])
    return np.dot(rot,pair)

def rotate_and_shift(pair,degrees,shift):
    return shift(_rotate(pair,degrees), shift)

def split(pair,degrees,shift):
    return rotate(_shift(_rotate(pair,degrees), shift),-degrees)

def unsplit(pair,degrees,shift):
    return split(pair,degrees+90,shift)
    
def taper(pair,width,centre=None):
    """Apply Hanning window about c0 sample
       of width number of samples, truncates traces"""
    
    if centre is None:
        centre = math.floor(pair.shape[1]/2)
        
    if width > pair.shape[1]:
        raise Exception('taper width is greater than trace length')
        
    t0 = centre - math.floor(width/2)
    t1 = centre + math.ceil(width/2)
    
    if t0 < 0:
        raise Exception('window starts before trace data')
    elif t1 > pair.shape[1]:
        raise Exception('window ends after trace data')
        
    return np.hanning(width) * pair[:,t0:t1]

def eigcov(pair):
    return np.sort(np.linalg.eigvals(np.cov(pair,rowvar=True)))
#     return np.sort(np.linalg.eigvals(np.cov(pair)))

def grideigcov(pair,maxshift,window=None, stepang=None,stepshift=None):
    
    # set some defaults
    if stepshift is None:
        stepshift = int(np.max([1,maxshift/40]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        window = int(np.min([pair.shape[1] * 0.5,maxshift * 10]))
    
    deg, shift = np.meshgrid(np.arange(0,180,stepang),
                             np.arange(0,maxshift,stepshift).astype(int))
    
    shape = deg.shape
    lam1 = np.zeros(shape)
    lam2 = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = rotate(pair,deg[0,ii]+90)
        for jj in np.arange(shape[0]):
            temp2 = shift(temp,shift[jj,ii])
            temp3 = taper(temp2,window)
            lam2[jj,ii], lam1[jj,ii] = eigcov(temp3)
    return deg, shift, lam1, lam2

def xcorr(pair):
    return np.correlate(pair[0,:],pair[1,:])[0]

def gridxcorr(pair,maxshift,window=None, stepang=None,stepshift=None):
    
    # set some defaults
    if stepshift is None:
        stepshift = int(np.max([1,maxshift/40]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        window = int(np.min([pair.shape[1] * 0.5,maxshift * 10]))
        
    deg, shift = np.meshgrid(np.arange(0,180,stepang),
                             np.arange(0,maxshift,stepshift).astype(int))
    
    shape = deg.shape
    xc = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = rotate(pair,deg[0,ii]+90)
        for jj in np.arange(shape[0]):
            temp2 = shift(temp,shift[jj,ii])
            temp3 = taper(temp2,window)
            xc[jj,ii] = xcorr(temp3)
    return deg, shift, xc  
    