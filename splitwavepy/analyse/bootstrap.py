"""
Use Bootstrapping to estimate measurement errors following (Sandvol and Hearn, 1994)
"""

from .pair import Pair
from .core import core
from .eigen import EigenM

import numpy as np

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

def bootstrap_sample(data,fast,lag):    
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
    
def bootstrap_loop(data,N=50):
    """
    Return list of bootstrap samples
    """        
    # initial measurement:
    m = sw.EigenM(data)
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
    bslist = [ sw.EigenM(bs) for bs in [ bootstrap_sample(data,degs,lags) for lags,degs in zip(m.tlags[idx],m.degs[idx]) ] ]
    return bslist
    
def boot_std(listM):
    avg = np.average(np.stack(listM))
    diffs = [ M - avg for M in listM ]
    np.transpose(diffs)
    ( 1 / ( len(listM) - 1 ) ) * [ np.transpose(m - avg)]
    