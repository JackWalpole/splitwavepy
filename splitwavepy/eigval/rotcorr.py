# -*- coding: utf-8 -*-
"""
The eigenvalue method of Silver and Chan (1991)
Low level routines works on numpy arrays and shifts using samples (doesn't know about time)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
from ..core.window import Window

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, stats
import math


     
def gridxc(x, y, degs, slags, window, **kwargs):
    """
    Grid search for splitting parameters applied to data.
    
    x is srcpol
    y is transverse
    lags = 1-D array of sample shifts to search over, if None an attempt at finding sensible values is made
    degs = 1-D array of rotations to search over, if None an attempt at finding sensible values is made
    window = Window object (if None will guess an appropriate window)
    rcvcorr = receiver correction parameters in tuple (fast,lag) 
    srccorr = source correction parameters in tuple (fast,lag) 
    """
    
                                 
    # grid of degs and lags to search over
    degs, lags = np.meshgrid(degs,slags)
    shape = degs.shape
    xc = np.zeros(shape)
    
    # avoid using "dots" in loops for performance
    rotate = core.rotate
    lag = core.lag
    chop = core.chop
    
    # pre-apply receiver correction
    if 'rcvcorr' in kwargs:
        x,y = core.unsplit(x,y,*kwargs['rcvcorr'])
    
    # make function to do source correction (used in loop)
    if 'srccorr' in kwargs:
        srcphi, srclag = kwargs['srccorr']
        def srccorr(x,y,ang):
            # unwind rotation
            x,y = rotate(x,y,srcphi-ang)
            # remove splitting
            x,y = lag(x,y,-srclag)
            return x,y
    else:
        def srccorr(x,y,ang):
            # no source correction so do nothing
            return x,y
    
    for ii in np.arange(shape[1]):
        tx, ty = rotate(x,y,degs[0,ii])
        for jj in np.arange(shape[0]):
            # remove splitting so use inverse operator (negative lag)
            ux, uy = lag(tx,ty,-lags[jj,ii])
            # if requested -- post-apply source correction
            ux, uy = srccorr(ux,uy,degs[0,ii])
            # chop to analysis window
            ux, uy = chop(ux,uy,window=window)
            # measure cross-correlation coefficient
            norm = np.sqrt(np.sum(ux**2) * np.sum(uy**2))
            xc[jj,ii] = np.correlate(ux,uy)/norm
            
    return degs,lags,xc


def ndf(y,window=None,detrend=False):
    """
    Estimates number of degrees of freedom using noise trace y.
    Uses the improvement found by Walsh et al (2013).
    """
        
    if detrend is True:
        # ensure no trend on the noise trace
        y = signal.detrend(y)

    if window is not None:
        # chop trace to window limits
        y = core.chop(y,window=window)
  
    Y = np.fft.fft(y)
    amp = np.absolute(Y)
    
    # estimate E2 and E4 following Walsh et al (2013) 
    a = np.ones(Y.size)
    a[0] = a[-1] = 0.5
    
    # equation (25)
    E2 = np.sum( a * amp**2)
    # equation (26)
    E4 = np.sum( (4 * a**2 / 3) * amp**4)
    
    # equation (31)
    ndf = 2 * ( 2 * E2**2 / E4 - 1 )
    
    return ndf
    
def ftest(xc,ndf,alpha=0.05):
    """
    returns (positive) xc value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Wuestefeld et al. (2008)
    """
    
    # ensure xc values are positive
    xcmax = np.max(np.abs(xc))
    
    # Fisher Transform
    z = math.atanh(xcmax)
    
    # Lookup F value
    k = 2 # two parameters, phi and dt.
    F = stats.f.ppf(1-alpha,k,ndf)
    
    # 
    # z_crit = z - ((z * k ) / ((ndf-k) * F)) * math.sqrt(1/(ndf-3))
    z_crit = z * ( 1 + (k/(ndf-k)) * F ) * math.sqrt( 1 / (ndf-3) ) 
    xc_alpha = math.tanh(z_crit) 
    return xc_alpha
    
def conf95(xc,ndf):
    """
    returns (positive) xc value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Wuestefeld et al. (2008)
    """
    
    # ensure xc values are positive
    xcmax = np.max(np.abs(xc))
    
    # Fisher Transform
    z = math.atanh(xcmax)
    std = math.sqrt( 1 / (ndf-3) )
    z95 = z - 1.96 * std
    xc95 = math.tanh(z95) 
    return xc95
    
     
