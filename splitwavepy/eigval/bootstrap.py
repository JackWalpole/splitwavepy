"""
Use Bootstrapping to estimate measurement errors following (Sandvol and Hearn, 1994).

The method is generalised to work without prior knowlege of the wave source polarisation
and therefore does not can be used on phases other than SKS.
"""

from ..core.pair import Pair
from ..core import core
from .eigenM import EigenM

import numpy as np

class Bootstrap:
    
    def __init__(self,pair,**kwargs):
        self.listM = bootstrap_measurements(pair,**kwargs)
        v = [ m.lam1 / m.lam2 for m in self.listM ]
        self.stk_l1_l2 = np.stack(v)



def get_noise(y):
    """
    Return a randomly simulated noise trace with similar spectral properties to y.
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

def bootstrap_sample(data,fast,lag,**kwargs):
    """
    Return data with new noise sequence
    """    
    # copy original data
    bs = data.copy()   
    origang = bs.angle
    # replace noise sequence
    bs.unsplit(fast,lag)
    bs.rotateto(bs.pca())
    bs.y = get_noise(bs.y)
    bs.rotateto(origang)
    bs.split(fast,lag)
    return bs

def rho(n,step):
    """
    Polar density of measurements
    """
    if n == 0:
        return 1 / (np.pi/4 * step**2)
    elif n > 0:
        return 1 / (2 * np.pi * n * step)
    else:
        raise Exception('n not valid')
        
rho = np.vectorize(rho)
    
def bootstrap_measurements(data,N=50,**kwargs):
    """
    Return list of bootstrap measurements
    """        
    # initial measurement:
    m = EigenM(data,**kwargs)
    mlags = m.tlags[:,0]
    mdegs = m.degs[0,:]
    # get probability surface to pick from
    # boost surf by **3 to enhance probability of picks at peaks (value chosen by testing on synthetics)
    surf = (m.lam1/m.lam2)**3
    dlag = m.lags[1,0] - m.lags[0,0]
    density = rho(m.lags,dlag)
    surf = surf / density
    surf = surf / surf.sum()
    
    # pick fast and tlag from surf
    probs = surf.ravel()
    picks = np.random.choice(probs.size,size=N,replace=True,p=probs)
    idx = np.unravel_index(picks,surf.shape)
    
    # generate bootstrap sample measurements    
    bslist = [ EigenM(bs,tlags=mlags,degs=mdegs) for bs in [ bootstrap_sample(data,degs,lags) for lags,degs in zip(m.tlags[idx],m.degs[idx]) ] ]
    return bslist
    
# def boot_std(listM):
#     avg = np.average(np.stack(listM))
#     diffs = [ M - avg for M in listM ]
#     np.transpose(diffs)
#     ( 1 / ( len(listM) - 1 ) ) * [ np.transpose(m - avg)]
    