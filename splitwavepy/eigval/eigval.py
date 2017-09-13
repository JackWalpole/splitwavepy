"""
The eigenvalue method of Silver and Chan (1991)
Low level routines works on numpy arrays and shifts using samples (doesn't know about time)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# import core as c
# from . import core as c
# from . import plotting as p

import .core

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, stats


    
def eigvalcov(data):
    """
    return sorted eigenvalues of covariance matrix
    lambda1 first, lambda2 second
    """
    return np.sort(np.linalg.eigvals(np.cov(data)))
    
    
def grideigval(data, maxlag=None, window=None, stepang=None, steplag=None):

    # set some defaults
    if maxlag is None:
        maxlag = int(data[0].size / 10)
        maxlag = maxlag if maxlag%2==0 else maxlag + 1
    if steplag is None:
        steplag = 2 * int(np.max([1,maxlag/80]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        # ensure window is odd length
        window = int(np.min([data.shape[1] * 0.5,maxlag * 10]))
        window = window if window%2==1 else window + 1

    degs, lags = np.meshgrid(np.arange(0,180,stepang),
                             np.arange(0,maxlag,steplag).astype(int))

    shape = degs.shape
    lam1 = np.zeros(shape)
    lam2 = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = core.rotate(data,degs[0,ii])
        for jj in np.arange(shape[0]):
            # remove splitting so use inverse operator (negative lag)
            temp2 = core.lag(temp,-lags[jj,ii])
            temp3 = core.window(temp2,window)
            lam2[jj,ii], lam1[jj,ii] = eigvalcov(temp3)
            
    return data,degs,lags,lam1,lam2,window

# def grideigval(data, maxlag=None, window=None, stepang=None, steplag=None):
#     data,degs,lags,lam1,lam2,window = _grideigval(data)
#     return Measurement(data=data,degs=degs,lags=lags,lam1=lam1,lam2=lam2,window=window)
    

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
    
def ftest(lam2min,ndf,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Silver and Chan (1991)
    """
    # lam2min = lam2.min()
    k = 2 # two parameters, phi and dt.
    # R = ((lam2 - lam2min)/k) /  (lam2min/(ndf-k))
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha
    

     
