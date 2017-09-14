"""
The eigenvalue method of Silver and Chan (1991)
Low level routines works on numpy arrays and shifts using samples (doesn't know about time)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
# from . import plotting

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, stats


    
def eigvalcov(data):
    """
    return sorted eigenvalues of covariance matrix
    lambda1 first, lambda2 second
    """
    return np.sort(np.linalg.eigvals(np.cov(data)))
    
    
def grideigval(data, lags=None, degs=None, window=None,rcvcorr=None,srccorr=None):
    """
    Grid search for splitting parameters applied to data.
    
    lags = 1-D array of sample shifts to search over, if None an attempt at finding sensible values is made
    degs = 1-D array of rotations to search over, if None an attempt at finding sensible values is made
    window = number of samples included in analysis (must be odd)
    rcvcorr = receiver correction parameters in tuple (fast,lag) 
    srccorr = source correction parameters in tuple (fast,lag) 
    """
    # set some defaults
    if lags is None:
        maxlag = int(data[0].size / 10)
        maxlag = maxlag if maxlag%2==0 else maxlag + 1
        steplag = 2 * int(np.max([1,maxlag/80]))
        lags = np.arange(0,maxlag,steplag).astype(int)        
        
    if degs is None:
        stepang = 2
        degs = np.arange(0,180,stepang)
        
    if window is None:
        # by default whatevers smaller: half trace length or 10 * maxlag
        # ensure window is odd length
        window = int(np.min([data.shape[1] * 0.5,maxlag * 10]))
        window = window if window%2==1 else window + 1

    # if requested -- pre-apply receiver correction
    if rcvcorr is not None:
        data = core.unsplit(data,*rcvcorr)
        
    gdegs, glags = np.meshgrid(degs,lags)

    shape = gdegs.shape
    lam1 = np.zeros(shape)
    lam2 = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = core.rotate(data,gdegs[0,ii])
        for jj in np.arange(shape[0]):
            # remove splitting so use inverse operator (negative lag)
            temp2 = core.lag(temp,-glags[jj,ii])
            # if requested -- post-apply source correction
            if srccorr is not None:
                temp2 = core.unsplit(temp2,*srccorr)
            temp2 = core.window(temp2,window)
            lam2[jj,ii], lam1[jj,ii] = eigvalcov(temp2)
            
    return gdegs,glags,lam1,lam2,window

def ndf(y,taper=False,detrend=True):
    """
    Estimates number of degrees of freedom using noise trace y.
    Uses the improvement found by Walsh et al (2013).
    By default will detrend data to ensure zero mean
    and will taper edges using a Tukey filter affecting amplitudes of data at edges (extreme 5%)
    """

    if taper is True:
        y = y * signal.tukey(y.size,0.05)
        
    if detrend is True:
        # ensure no trend on the noise trace
        y = signal.detrend(y)

  
    Y = np.fft.fft(y)
    amp = np.absolute(Y)
    
    # estimate E2 and E4 following Walsh et al (2013)
    a = np.ones(Y.size)
    a[0] = a[-1] = 0.5
    E2 = np.sum( a * amp**2)
    E4 = (np.sum( (4 * a**2 / 3) * amp**4))
    
    ndf = 2 * ( 2 * E2**2 / E4 - 1 )
    
    return ndf
    
def ftest(lam2,ndf,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Silver and Chan (1991)
    """
    lam2min = lam2.min()
    k = 2 # two parameters, phi and dt.
    # R = ((lam2 - lam2min)/k) /  (lam2min/(ndf-k))
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha
    

     