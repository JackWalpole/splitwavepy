"""
The eigenvalue method of Silver and Chan (1991)
"""

from . import core as c
from . import plotting as p
import numpy as np
import matplotlib.pyplot as plt

class Measurement:
    
    def __init__(self,degs=None,lags=None,lam1=None,lam2=None):
        
        self.method = 'Eigenvalue'
        
        if degs is None or lags is None or lam1 is None or lam2 is None:
            # generate synthetic
            self.degs, self.lags, self.lam1, self.lam2 = _synthM()     
        else:
            self.degs = degs
            self.lags = lags
            self.lam1 = lam1
            self.lam2 = lam2
            

        
    def min_idx(self,vals):
        """
        return indices of min value in vals grid
        """
        return np.unravel_index(np.argmin(vals),vals.shape)
        
    def max_idx(self,vals):
        """
        return indice of max value in vals grid
        """
        return np.unravel_index(np.argmax(vals),vals.shape)
    

    def plot(self,vals=None,cmap='magma'):
        """
        plot the measurement.
        by default plots lam1/lam2
        """
        if vals is None:
            vals = self.lam1 / self.lam2
            
        p.plot_surf(self.lags,self.degs,vals,cmap=cmap)        

        
        
    
def eigvalcov(data):
    """
    return sorted eigenvalues of covariance matrix
    lambda1 first, lambda2 second
    """
    return np.sort(np.linalg.eigvals(np.cov(data)))
    
def eigcov(data):
    """
    return eigenvalues and eigenvectors of covariance matrix
    """
    vals,vecs = np.linalg.eig(np.cov(data))
    
    
def _grideigvalcov(data, maxshift=None, window=None, stepang=None, stepshift=None):

    # set some defaults
    if maxshift is None:
        maxshift = int(data[0].size / 10)
        maxshift = maxshift if maxshift%2==0 else maxshift + 1
    if stepshift is None:
        stepshift = 2 * int(np.max([1,maxshift/20]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        # ensure window is odd length
        window = int(np.min([data.shape[1] * 0.5,maxshift * 10]))
        window = window if window%2==1 else window + 1

    degs, lags = np.meshgrid(np.arange(0,180,stepang),
                             np.arange(0,maxshift,stepshift).astype(int))

    shape = degs.shape
    lam1 = np.zeros(shape)
    lam2 = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = c.rotate(data,degs[0,ii])
        for jj in np.arange(shape[0]):
            # remove splitting so use inverse operator (negative lag)
            temp2 = c.lag(temp,-lags[jj,ii])
            temp3 = c.window(temp2,window)
            lam2[jj,ii], lam1[jj,ii] = eigvalcov(temp3)
            
    return degs,lags,lam1,lam2

def grideigvalcov(data, maxshift=None, window=None, stepang=None, stepshift=None):
    degs,lags,lam1,lam2 = _grideigvalcov()
    return Measurement(degs=degs,lags=lags,lam1=lam1,lam2=lam2)

def ndf(y,taper=True,detrend=True):
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
    
def Ftest(M,ndf,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha=0.05 yielding the value of lambda2 at 95% confidence interval
    """
    lam2min = M.lam2.min()
    k = 2
    R = ((M.lam2 - lam2min)/k) /  (lam2min/(ndf-k))
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha
    
def _synthM(deg=25,lag=10):
    P = c.Pair()
    P.split(deg,lag)
    return _grideigvalcov(P.data)
     
