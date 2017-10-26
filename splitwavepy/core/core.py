# -*- coding: utf-8 -*-
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
import math

##############

def near(x): return np.rint(x).astype(int)
def even(x): return 2*np.rint(x/2).astype(int)    
def odd(x): return (2*np.rint(np.ceil(x/2))-1).astype(int)

def time2samps(t,delta,mode='near'):
    """
    convert a time to number of samples given the sampling interval.
    """
    rat = (t / delta)  
    if mode == 'near': return near(rat)
    if mode == 'even': return even(rat)
    if mode == 'odd' : return odd(rat)

def samps2time(samps,delta):
    """
    convert a number of samples to time given the sampling interval.
    """
    return samps * delta
    
################

def lag(x,y,samps):
    """
    Lag x samps/2 to the left and
    lag y samps/2 to the right.
    samps must be even.
    If samps is negative x shifted to the right and y to the left.
    This process truncates trace length by samps and preserves centrality.
    Therefore windowing must be used after this process 
    to ensure even trace lengths when measuring splitting.
    """
    if samps == 0:
        return x,y

    if samps > 0:
        # positive shift
        return x[samps:], y[:-samps]
    else:
        # negative shift
        return x[:samps], y[-samps:]
      
def rotate(x,y,degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates from x to y axis
       e.g. N to E if row 0 is N cmp and row1 is E cmp"""
    ang = math.radians(degrees)
    rot = np.array([[ np.cos(ang), np.sin(ang)],
                    [-np.sin(ang), np.cos(ang)]])
    xy = np.dot(rot, np.vstack((x,y)))
    return xy[0], xy[1]

def split(x,y,degrees,samps):
    """Apply forward splitting and rotate back"""
    if samps == 0:
        return x,y
    x,y = rotate(x,y,degrees)
    x,y = lag(x,y,samps)
    x,y = rotate(x,y,-degrees)
    return x,y

def unsplit(x,y,degrees,samps):
    """Apply inverse splitting and rotate back"""
    return split(x,y,degrees,-samps)

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
   
def eigcov(data):
    """
    Return eigen values and vectors of covariance matrix
    """
    eigenValues, eigenVectors = np.linalg.eig(np.cov(data))
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return eigenValues, eigenVectors
    
    
# def snr(data):
#     """
#     Returns signal to noise ratio assuming signal on trace1 and noise on trace2
#     Calculates on each trace by sum of squares and then takes the ratio
#     """samps
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
    """Gaussian noise convolved with a (normalised) gaussian wavelet.
       samps = size,
       sigma  = amp,
       width of gaussian = smooth.
    """
    norm = 1/(smooth*np.sqrt(2*np.pi))
    gauss = norm * signal.gaussian(size,smooth)
    n = np.random.normal(0,amp,size)
    return np.convolve(n,gauss,'same')  
    
def resample_noise(y):
    """
    Return a randomly simulated noise trace with similar spectral properties to y.
    
    Following Sandvol and Hearn.
    """  
    # white noise
    x = np.random.normal(0,1,y.size)
    # convolve with y
    x = np.convolve(x,y,'same')
    # additional randomisation
    x = np.roll(x,np.random.randint(y.size))
    # whipeout near nyquist
    x = np.convolve(np.array([1,1,1]),x,'same')
    # normalise energy
    x = x * np.sqrt((np.sum(y**2) / np.sum(x**2)))
    # return
    return x
    
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
