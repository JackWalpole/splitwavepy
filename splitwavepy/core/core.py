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
from scipy import signal, stats
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
    Lag x samps to the left and
    lag y samps to the right.
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

## Measurement 
   
def eigcov(data):
    """
    Return eigen values and vectors of covariance matrix
    """
    eigenValues, eigenVectors = np.linalg.eig(np.cov(data))
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return eigenValues, eigenVectors
    
# def eigvalcov(data):
#     """
#     return sorted eigenvalues of covariance matrix
#     lambda2 first, lambda1 second
#     """
#     return np.sort(np.linalg.eigvalsh(np.cov(data)))
    
def eigvalcov(x,y):
    """
    return sorted eigenvalues of covariance matrix
    lambda2 first, lambda1 second
    """
    data = np.vstack((x,y))
    return np.sort(np.linalg.eigvalsh(np.cov(data)))
  
def transenergy(x,y):
    """
    return energy
    lambda1 first, lambda2 second
    """
    energy = lambda x: np.sum(x**2)
    return energy(x), energy(y) 
    
def crosscorr(x,y):
    norm = math.sqrt(np.sum(x**2) * np.sum(y**2))
    xc = np.correlate(x,y)/norm
    return xc

def crossconv(obsx, obsy, prex, prey):
    """
    Cross convolve 
    """
    x = np.convolve(obsx, prey)
    y = np.convolve(prex, obsy)
    return x, y

def misfit(x, y):
    num = np.trapz((x - y)**2)
    den = np.trapz(x**2) + np.trapz(y**2)
    return num / den  
    
def crossconvmf(obsx, obsy, prex, prey):
    x, y = crossconv(obsx, obsy, prex, prey)
    return misfit(x, y)  

def splittingintensity(rad,trans):
    """
    Calculate splitting intensity.
    """    
    rdiff = np.gradient(rad)
    s = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
    return s

# Errors

def ndf(y):
    """
    Estimates number of degrees of freedom using noise trace y.
    Uses the improvement found by Walsh et al (2013).
    """
  
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
    
def ftest(lam2,ndf,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Silver and Chan (1991)
    """
    
    # check ndf is big enough
    if ndf < 3:
        raise Exception('Number of degrees of freedom is less than 3.  This likely indicates a problem which would lead to a spurios mesaurement.  Check window length.')
    
    lam2min = lam2.min()
    k = 2 # two parameters, phi and dt.
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha    

# Null Criterion

def Q(fastev,lagev,fastrc,lagrc):
    """Following Wuestefeld et al. 2010"""
    omega = math.fabs((fastev - fastrc + 3645)%90 - 45) / 45
    delta = lagrc / lagev
    dnull = math.sqrt(delta**2 + (omega-1)**2) * math.sqrt(2)
    dgood = math.sqrt((delta-1)**2 + omega**2) * math.sqrt(2)
    if dnull < dgood:
        return -(1 - dnull)
    else:
        return (1 - dgood)

# Signal to noise

def snrRH(x,y):
    """
    Returns signal to noise ratio assuming signal on trace1 and noise on trace2
    Uses the method of Restivo and Helffrich (1999):
    peak amplitude on trace1 / 2*std trace2
    """
    signal = np.max(x)
    noise = 2 * np.std(y)
    return signal / noise

# Useful bits and pieces

def synth(**kwargs):
    """return ricker wavelet synthetic data"""
    
    # defaults    
    if 'pol' not in kwargs: kwargs['pol'] = 0.
    if 'delta' not in kwargs: kwargs['delta'] = 1.
    if 'split' not in kwargs: kwargs['split'] = []
    if 'noise' not in kwargs: kwargs['noise'] = 0.001
    if 'nsamps' not in kwargs: kwargs['nsamps'] = 1001
    if 'width' not in kwargs: kwargs['width'] = 32
    if 'noisewidth' not in kwargs: kwargs['noisewidth'] = kwargs['width']/4

    # initiate data with clean ricker wavelet
    nsamps = int(kwargs['nsamps'])  
    x = signal.ricker(nsamps, kwargs['width'])
    y = np.zeros(nsamps)
    
    # rotate to polarisation 
    # negative because we are doing the active rotation of data, whereas
    # core.rotate does the passive transormation of the co-ordinate system
    x,y = rotate(x, y, -kwargs['pol'])

    if isinstance(kwargs['split'], tuple):
        fast, lag = kwargs['split']
        # add any splitting -- lag samples must be even
        slag = time2samps(lag, kwargs['delta'], mode='even')
        x,y = split(x, y, fast, slag)
    elif isinstance(kwargs['split'], list):        
        for parms in kwargs['split']:
            fast, lag = parms
            # add any splitting -- lag samples must be even
            slag = time2samps(lag, kwargs['delta'], mode='even')
            x,y = split(x, y, fast, slag)
    
    # add noise - do this last to avoid splitting the noise
    x = x + noise(x.size, kwargs['noise'], int(kwargs['noisewidth']))    
    y = y + noise(x.size, kwargs['noise'], int(kwargs['noisewidth']))

    return x,y
    
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
