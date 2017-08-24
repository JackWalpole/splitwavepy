"""
The eigenvalue method of Silver and Chan (1991)
"""

from . import core as c
from . import plotting as p
import numpy as np
import matplotlib.pyplot as plt

class Measurement:
    
    def __init__(self,degs,lags,lam1,lam2):
        self.method = 'Eigenvalue'
        self.degs = degs
        self.lags = lags
        self.lam1 = lam1
        self.lam2 = lam2
        self.phi = np.argmin()
        self.lag = np.argmin()        
        
    # def print_summary():
    #     """
    #     Print useful info about the measurement
    #     """
    #     print()
    
    # methods
    def plot_val(self,vals=self.lam1/self.lam2):
        """
        plot a surface
        """
        plt.contourf(lags,degs,vals,cmap='viridis')
        plt.show()
        
    def min_idx(self,vals):
        return np.unravel_index(np.argmin(vals),vals.shape)
        
    def max_idx(self,vals):
        return np.unravel_index(np.argmax(vals),vals.shape)
    

        
        
    
def eigvalcov(pair):
    """
    return sorted eigenvalues of covariance matrix
    lambda1 first, lambda2 second
    """
    return np.sort(np.linalg.eigvals(np.cov(pair)))
    
def eigcov(pair):
    """
    return eigenvalues and eigenvectors of covariance matrix
    """
    vals,vecs = np.linalg.eig(np.cov(pair))
    
    
def grideigvalcov(pair, maxshift=None, window=None, stepang=None, stepshift=None):

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
            lam2[jj,ii], lam1[jj,ii] = eigvalcov(temp3)
            
    return Measurement(deg,lag,lam1,lam2)




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
    
    
