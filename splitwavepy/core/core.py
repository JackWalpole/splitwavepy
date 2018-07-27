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

# from .window import Window

import numpy as np
from scipy import signal, stats
from scipy.interpolate import interp1d 
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

def lag(x, y, samps):
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
        return x, y

    if samps > 0:
        # positive shift
        return x[samps:], y[:-samps]
    else:
        # negative shift
        return x[:samps], y[-samps:]
      
# def rotate(x,y,degrees):
#     """row 0 is x-axis and row 1 is y-axis,
#        rotates from x to y axis
#        e.g. N to E if row 0 is N cmp and row1 is E cmp"""
#     ang = math.radians(degrees)
#     rot = np.array([[ np.cos(ang), np.sin(ang)],
#                     [-np.sin(ang), np.cos(ang)]])
#     xy = np.dot(rot, np.vstack((x,y)))
#     return xy[0], xy[1]

def rotate(x, y, degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates from x to y axis
       e.g. N to E if row 0 is N cmp and row1 is E cmp"""
    ang = math.radians(degrees)
    cang = math.cos(ang)
    sang = math.sin(ang)
    rot = np.array([[ cang, sang],
                    [-sang, cang]])
    xy = np.dot(rot, np.vstack((x, y)))
    return xy[0], xy[1]
    
def _rot(degrees):
    ang = math.radians(degrees)
    cang = math.cos(ang)
    sang = math.sin(ang)
    rot = np.array([[ cang, sang],
                    [-sang, cang]])
    return rot
    
# def rotrad(x, y, radians):
#     """row 0 is x-axis and row 1 is y-axis,
#        rotates from x to y axis
#        e.g. N to E if row 0 is N cmp and row1 is E cmp"""
#     cang = math.cos(radians)
#     sang = math.sin(radians)
#     rot = np.array([[ cang, sang],
#                     [-sang, cang]])
#     xy = np.dot(rot, np.vstack((x, y)))
#     return xy[0], xy[1]

def split(x, y, degrees, samps):
    """Apply forward splitting and rotate back"""
    if samps == 0:
        return x, y
    x, y = rotate(x, y, degrees)
    x, y = lag(x, y, samps)
    x, y = rotate(x, y, -degrees)
    return x, y

def unsplit(x, y, degrees, samps):
    """Apply inverse splitting and rotate back"""
    return split(x, y, degrees, -samps)
    
def chop(x, s0, s1):
    """Chop 1-d numpy arrays from s0 to s1"""
    return x[s0:s1]
    
def taper(x, alpha=0.9):
    """Taper data in x using a Tukey window, also known as a tapered cosine window."""
    return x * signal.tukey(x.size, alpha=alpha)


# def chop(*args,**kwargs):
#     """Chop trace, or traces, using window"""
#
#     if ('window' in kwargs):
#         window = kwargs['window']
#
#     if not isinstance(window,Window):
#         raise Exception('window must be a Window')
#
#     length = args[0].size
#
#     if window.width > length:
#         raise Exception('window width is greater than trace length')
#
#     centre = int(length/2) + window.offset
#     hw = int(window.width/2)
#     t0 = centre - hw
#     t1 = centre + hw
#
#     if t0 < 0:
#         raise Exception('chop starts before trace data')
#     elif t1 > length:
#         raise Exception('chop ends after trace data')
#
#     if window.tukey is not None:
#         tukey = signal.tukey(window.width,alpha=window.tukey)
#     else:
#         tukey = 1.
#
#     if len(args)==1:
#         return args[0][t0:t1+1] * tukey
#     elif len(args)==2:
#         return args[0][t0:t1+1] * tukey, args[1][t0:t1+1] * tukey
#     elif len(args)==3:
#         return args[0][t0:t1+1] * tukey, args[1][t0:t1+1] * tukey, args[2][t0:t1+1] * tukey

## Measurement 
   
def eigcov(x, y):
    """
    Return eigen values and vectors of covariance matrix
    """
    data = np.vstack((x, y))
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
    
def eigvalcov(x, y, **kwargs):
    """
    return sorted eigenvalues of covariance matrix
    lambda1 first, lambda2 second
    """
    data = np.vstack((x,y))
    return np.sort(np.linalg.eigvalsh(np.cov(data)))[::-1]
  
def transenergy(x, y, **kwargs):
    """
    return energy
    lambda1 first, lambda2 second
    """
    # energy = lambda x: np.sum(x**2)
    energy = lambda x: np.var(x)
    return energy(x), energy(y) 
    
def crosscorr(x, y, **kwargs):
    norm = math.sqrt(np.sum(x**2) * np.sum(y**2))
    xc = np.correlate(x, y)/norm
    return xc
    
def pearson(x, y, **kwargs):
    x = x - np.mean(x)
    y = y - np.mean(y)
    norm = math.sqrt(np.sum(x**2) * np.sum(y**2))
    xc = np.correlate(x, y)/norm
    return xc
    
def fisher_r(rho):
    return np.arctanh(rho)

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

# def splittingintensity(rad, trans):
#     """
#     Calculate splitting intensity.
#     """
#     rdiff = np.gradient(rad)
#     s = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
#     return s

# Grid covariance

def slagchop(x, y, w0, w1, slag):
    d = int(slag/2)
    return x[w0+d: w1+d], y[w0-d: w1-d] 

def running_mean(x, w0, w1, slags):
    d = int(slags[-1]/2)
    x = x[w0-d: w1+d]
    n = w1-w0
    return np.convolve(x, np.ones((n,))/n, mode='valid')
    
# def gridcov(x, y, w0, w1, degs, slags):
#     # prepare a list of data rotated to degs
#     rot_data = [ rotate(x, y, deg) for deg in degs ]
#     # prepare empty covariance arrays
#     gridcov = np.empty((degs.size, slags.size, 2, 2))
#     c = np.empty((2, 2))
#     ii = 0
#     # now loop and calculate
#     for rot in rot_data:
#         # this is the mean in each window
#         meanx = running_mean(rot[0], w0, w1, slags)
#         meany = running_mean(rot[1], w0, w1, slags)
#         jj = 0
#         for slag in slags:
#             wx, wy  = slagchop(*rot, w0, w1, slag)
#             dx, dy = wx - meanx[slag], wy - meany[slag]
#             n = dx.size
#             c[0, 0] = np.sum(dx * dx)
#             c[1, 0] = c[0, 1] = np.sum(dx * dy)
#             c[1, 1] = np.sum(dy * dy)
#             c = c / n
#             gridcov[ii, jj, :, :] = c
#             jj += 1
#         ii += 1
#     return gridcov

# def gridcov(x, y, w0, w1, degs, slags):
#     # prepare empty covariance arrays
#     g = np.empty((degs.size, slags.size, 2, 2))
#     n = w1 - w0
#     npsum = np.sum # remove dots from inner loop
#     # now loop and calculate
#     for ii in range(degs.size):
#         # prepare a list of data rotated to degs
#         rot = rotate(x, y, degs[ii])
#         # this is the mean in each window
#         meanx = running_mean(rot[0], w0, w1, slags)
#         meany = running_mean(rot[1], w0, w1, slags)
#         # loop over lags
#         for jj in range(slags.size):
#             slag = slags[jj]
#             wx, wy  = slagchop(rot[0], rot[1], w0, w1, -slag)
#             dx, dy = wx - meanx[slag], wy - meany[slag]
#             g[ii, jj, 0, 0] = npsum(dx * dx)
#             g[ii, jj, 1, 0] = g[ii, jj, 0, 1] = npsum(dx * dy)
#             g[ii, jj, 1, 1] = npsum(dy * dy)
#     return g / n

def gridcov(x, y, w0, w1, degs, slags):
    # prepare empty covariance arrays
    g = np.empty((slags.size, degs.size, 2, 2))
    n = w1 - w0
    npsum = np.sum # remove dots from inner loop
    # now loop and calculate
    for ii in range(degs.size):
        # prepare a list of data rotated to degs
        rot = rotate(x, y, degs[ii])
        # this is the mean in each window
        meanx = running_mean(rot[0], w0, w1, slags)
        meany = running_mean(rot[1], w0, w1, slags)
        # loop over lags
        for jj in range(slags.size):
            slag = slags[jj]
            wx, wy  = slagchop(rot[0], rot[1], w0, w1, -slag)
            dx, dy = wx - meanx[slag], wy - meany[slag]
            g[jj, ii, 0, 0] = npsum(dx * dx)
            g[jj, ii, 1, 0] = g[jj, ii, 0, 1] = npsum(dx * dy)
            g[jj, ii, 1, 1] = npsum(dy * dy)
    return g / n    



def gridcovfreq(x, y, ndegs=90, nslags=50):
    """Returns grid of covariance matrices of shape (ndegs/2, (2*nslags)-1, 2, 2).
       Use covfreq_reshape to get this into a more user friendly shape.
       x and y are the windowed traces.
       Probably best to use a reasonably long, tapered, window as shifts
       are performed implicitly in the Frequency domain (so wrap around?).
    """
    mdegs = int(ndegs/2)
    mlags = int(nslags*2) - 1
    degs = np.linspace(0, 90, mdegs, endpoint=False)
    n = x.size
    x = x - np.mean(x)
    y = y - np.mean(y)
    # cov 0 (calculate in time domain)
    sumxx = np.sum(x*x)
    sumyy = np.sum(y*y)
    sumxy = np.sum(x*y)
    ncov0 = [[sumxx, sumxy], [sumxy, sumyy]]
    g = np.empty((mlags, mdegs, 2, 2))
    # Fourier Transform
    fx = np.fft.fft(x)
    fy = np.fft.fft(y)
    fxy = np.vstack((fx, fy))
    # now loop and calculate
    for ii in range(mdegs):
        # rotate
        rot = _rot(degs[ii])
        fxyr = np.dot(rot, fxy)
        fxr = fxyr[0]
        fyr = fxyr[1]
        # correlate
        cxy = fxr * fyr.conj()
        # inverse transform
        icxy = np.fft.ifft(cxy).real
        # fill g (the grid of covariance matrices)
        g[:,ii] = np.dot(rot, np.dot(ncov0, rot.T))
        covxy = np.roll(icxy, nslags-1)[0:mlags]
        # fill in cross-covariance between x and y
        g[:,ii,0,1] = g[:,ii,1,0] = covxy
    return g / n    


def cov_reshape(cov):
    """Reshape a covariance map to match the standard."""
    shp = cov.shape
    mid = int((shp[0]-1)/2)
    pos = cov[:mid+1,:,:,:]
    pos = np.flip(pos, 0)
    neg = cov[mid:,:,:,:]
    return np.concatenate((pos, neg), axis=1)


# def cov_rotate(cov, deg):
#     """Rotate covariance matrices to deg."""
#     shp = cov.shape
#     outcov = np.empty(shp)
#     degs = np.linspace(0, 90, shp[1], endpoint=False)
#     for ii in range(shp[1]):
#         rot = _rot(deg - degs[ii])
#         outcov[:, ii] = np.matmul(rot, np.matmul(cov[:,ii], rot.T))
#     return outcov

def cov_rotate(cov, deg):    
    """Rotate covariance matrices to deg."""
    shp = cov.shape
    outcov = np.empty(shp)
    degs = np.linspace(0, 180, shp[1], endpoint=False)
    for ii in range(shp[1]):
        rot = _rot(deg - degs[ii])
        outcov[:, ii] = np.matmul(rot, np.matmul(cov[:,ii], rot.T))
    return outcov

# def slagchop_srccorr(x, y, w0, w1, slag, srcfast, srcslag):
#     x, y = rot2(x, y, srcfast)
#     x, y = lag(x, y, srcslag)
#     x, y = rot2(x, y, -srcfast)
#     d = int(slag/2) - int(srcslag/2)
#     return x[w0+d: w1+d], y[w0-d: w1-d]
#
# def gridcov_srcorr(x, y, w0, w1, degs, slags, srcfast, srcslag):
#     # prepare a list of data rotated to degs
#     rot_data = [ rot2(x0, y0, deg) for deg in degs ]
#     # prepare empty covariance arrays
#     gridcov = np.empty((degs.size, slags.size, 2, 2))
#     c = np.empty((2, 2))
#     ii = 0
#     # now loop and calculate
#     for rot in rot_data:
#         jj = 0
#         for slag in slags:
#             wx, wy  = slagchop_srccorr(*rot, w0, w1, slag, srcfast, srcslag)
#             mx, my = np.mean(wx), np.mean(wy)
#             dx, dy = wx - mx, wy - my
#             n = dx.size
#             c[0, 0] = np.sum(dx * dx)
#             c[1, 0] = c[0, 1] = np.sum(dx * dy)
#             c[1, 1] = np.sum(dy * dy)
#             c = c / n
#             gridcov[ii, jj, :, :] = c
#             jj += 1
#         ii += 1
#     return gridcov

def covmap2eig(cov):
    eigvals, eigvecs = np.linalg.eigh(cov[:, :])
    return eigvals, eigvecs

def covmap2eigvals(cov):
    eigvals = np.linalg.eigvalsh(cov[:, :])
    lam2 = eigvals[:, :, 0]
    lam1 = eigvals[:, :, 1]
    return lam1, lam2
    
def covmap2rho(cov):
    stdx = np.sqrt(cov[:, :, 0, 0])
    stdy = np.sqrt(cov[:, :, 1, 1])
    rho = cov[:, :, 0, 1] / (stdx * stdy)
    return rho
    

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
    
def ftest(lam2, ndf, alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Silver and Chan (1991)
    """
    
    # check ndf is big enough
    if ndf < 3:
        raise Exception('Number of degrees of freedom is less than 3.  \
        This likely indicates a problem which would lead to a spurious measurement.  \
        Check window length.  Check data are demeaned.  Check frequency content.')
    
    lam2min = lam2.min()
    k = 2 # two parameters, phi and dt.
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha    
    
def val_at_alpha(data, alpha):
    """ Find value of function at the alpha level """
    idx = np.argsort(data) 
    cum = np.cumsum(data[idx])
    tot = np.max(cum)
    get_x_at_cum = interp1d(cum, np.arange(cum.size))
    get_val_at_x = interp1d(np.arange(data.size), data[idx])
    xval = get_x_at_cum(tot*alpha)
    return get_val_at_x(xval)
    
# Bootstrapping
    
def resample_noise(y):
    """
    Return a randomly simulated noise trace with similar spectral properties to y.
    
    Following Sandvol and Hearn.
    """  
    # white noise
    x = np.random.normal(0, 1, y.size)
    # convolve with y
    x = np.convolve(x, y, 'same')
    # additional randomisation
    x = np.roll(x, np.random.randint(y.size))
    # whipeout near nyquist
    x = np.convolve(np.array([1,1,1]), x, 'same')
    # normalise energy
    x = x * np.sqrt((np.sum(y**2) / np.sum(x**2)))
    # return
    return x
    
def bootstrap_resamp(x, y):
    """Resample data for bootstrapping"""
    idx = np.random.choice(x.size, x.size)
    return x[idx], y[idx]    
    
def kde(vals):
    # instantiate and fit the KDE model
    return stats.gaussian_kde(vals)

# Null Criterion

def q(fastsc, lagsc, fastxc, lagxc):
    """Following Wuestefeld et al. 2010"""
    omega = math.fabs((fastsc - fastxc + 3645)%90 - 45) / 45
    delta = lagxc / lagsc
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
    
def noise(size, amp, smooth):
    """Gaussian noise convolved with a (normalised) gaussian wavelet.
       samps = size,
       sigma  = amp,
       width of gaussian = smooth.
    """
    norm = 1/(smooth*np.sqrt(2*np.pi))
    gauss = norm * signal.gaussian(size,smooth)
    n = np.random.normal(0,amp,size)
    return np.convolve(n,gauss,'same')  



def min_idx(vals):
    """
    return indices of min value in vals grid
    """
    return np.unravel_index(np.argmin(vals), vals.shape)

def max_idx(vals):
    """
    return indice of max value in vals grid
    """
    return np.unravel_index(np.argmax(vals), vals.shape)
