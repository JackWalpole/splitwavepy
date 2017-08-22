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
    
def grideigcov(pair,maxshift,window=None, stepang=None,stepshift=None):

    # set some defaults
    if stepshift is None:
        stepshift = 2 * int(np.max([1,maxshift/40]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        # ensure window is odd length
        window = int(np.min([pair.shape[1] * 0.5,maxshift * 10]))
        if window%2 == 0:
            window = window + 1

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
    return deg, lag, lam1, lam2