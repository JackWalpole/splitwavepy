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
    x = np.random.normal(y.size)
    # convolve with y
    x = np.convolve(y,x,'same')
    # additional randomisation
    x = np.roll(np.random.randint(y.size))
    # whipeout near nyquist
    x = np.covolve(np.array([1,1,1]))
    # normalise energy
    x = x / np.sum(y**2)
    # return
    return x

def bootstrap_sample(data,fast,lag):    
    # copy original data
    bs = data.copy()     
    # replace noise sequence
    copy.unsplit(fast,lag)
    copy.y = get_noise(copy.y)
    copy.split(fast,lag)
    return copy
    
def bootstrap_loop(data,N=50):
    """
    Return list of bootstrap samples
    """
    
    if not isinstance(data,Pair):
        raise Exception('data must be a Pair')
    
    # initial measurement:
    m = EigenM(data)
    # polar normalise lam1/lam2 to use as weights
    surf = polnorm * m.lam1 / m.lam2
    del m

    # choose correction to find noise sequence
    idx = np.random.choice(surf.ravel())
    fast, lag =
    
    # bslist = []
    #
    # # Bootstrap Loop
    # for ii np.range(N):
    #
    #     # Generate Bootstrap Sample
    #     bs = bootstrap_sample(data,fast,lag)
    #
    #     # Measure
    #     bm = EigenM(bs)
    #
    #     # Add to list
    #     bslist = bslist.append(bm)
        
    bslist = [ EigenM(x) for bs in [ bootstrap_sample(x) for x in range(N) ] ]
    return bslist
    
def boot_std(listM):
    avg = np.average(np.stack(listM))
    diffs = [ M - avg for M in listM ]
    np.transpose(diffs)
    ( 1 / ( len(listM) - 1 ) ) * [ np.transpose(m - avg)]
    