"""
The eigenvalue method of Silver and Chan (1991)
"""

from . import core as c
from . import plotting as p
import numpy as np

class Measurement:
    
    def __init__(self,degs,lags,lam1,lam2):
        self.method = 'Eigenvalue'
        self.degs = degs
        self.lags = lags
        self.lam1 = lam1
        self.lam2 = lam2
        
    # def summary():
    #     """
    #     Print useful info about the measurement
    #     """
    #     print()
    
    # methods
    def plot(self):
        self.plot = p.plot_surf(self.lags,self.degs,self.lam1,self.lam2)
        
        
    
def eigcov(pair):
    """get eigenvalues of covariance matrix"""
    return np.sort(np.linalg.eigvals(np.cov(pair,rowvar=True)))
    
def grideigcov(pair, maxshift=None, window=None, stepang=None, stepshift=None):

    # set some defaults
    if maxshift is None:
        maxshift = int(pair[0].size / 10)
        maxshift = maxshift if maxshift%2==0 else maxshift + 1
    if stepshift is None:
        stepshift = 2 * int(np.max([1,maxshift/20]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        # ensure window is odd length
        window = int(np.min([pair.shape[1] * 0.5,maxshift * 10]))
        window = window if window%2==1 else window + 1


    deg, lag = np.meshgrid(np.arange(0,180,stepang),
                             np.arange(0,maxshift,stepshift).astype(int))

    shape = deg.shape
    lam1 = np.zeros(shape)
    lam2 = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = c.rotate(pair,deg[0,ii])
        for jj in np.arange(shape[0]):
            # remove splitting so use inverse operator (negative lag)
            temp2 = c.lag(temp,-lag[jj,ii])
            temp3 = c.window(temp2,window)
            lam2[jj,ii], lam1[jj,ii] = eigcov(temp3)
            
    return Measurement(deg,lag,lam1,lam2)


def ndf(y,taper=True,detrend=True,):
    """
    Uses the improvement found by Walsh et al (2013).
    By default will detrend data and taper the edges.
    """

    if taper is True:
        # taper edges to reduce transform artefacts
        y = y * signal.tukey(y.size,0.05)
        
    if detrend is True:
        # ensure no trend on the noise trace
        y = signal.detrend(y)

  
    Y = np.fft.fft(y)
    amp = np.absolute(Y)
    
    # estimate E2 and E4 following Walsh et al (2013)
    # note we do not scale the first and last samples by half before summing
    # as is done by Walsh following Silver and Chan.
    # this is because in practice we use a discrete fourier transform, 
    # and not a continuous transform for which the former approach would be correct.
    E2 = np.sum(amp**2)
    E4 = 4/3 * np.sum(amp**4)
    
    ndf = 2 * ( 2 * E2**2 / E4 - 1 )
    
    return int(round(ndf))
    
