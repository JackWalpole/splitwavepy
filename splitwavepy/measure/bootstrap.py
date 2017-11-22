# -*- coding: utf-8 -*-
"""
Use Bootstrapping to estimate measurement errors following (Sandvol and Hearn, 1994).

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core.pair import Pair
from ..core import core
from .eigenM import EigenM

import numpy as np

class Bootstrap:
    """
    Produce *nits* bootstrap samples (measurement objects)
    by resimulating noise sequence and remeasuring.
    
    keywords
    
    nits      number of bootstrap samples to produce
    """
    
    def __init__(self,pair,n=50,**kwargs):
        if not isinstance(pair,Pair):
            raise TypeError('expecting a pair')
        kwargs['n'] = n
        self.data = pair
        self.listM = bs_loop(pair,**kwargs)
        # self.stk_l1_l2 = np.stack([ m.lam1 / m.lam2 for m in self.listM ])
        # self.stk_fastprofile = np.stack([m.fastprofile() for m in self.listM ])
        # self.stk_lagprofile = np.stack([m.lagprofile() for m in self.listM ])

def bs_loop(pair,**kwargs):
    """
    Return list of bootstrap measurements
    """        
    # initial measurement:
    m = EigenM(pair,**kwargs)
    mlags = m.lags[:,0]
    mdegs = m.degs[0,:]
    # get probability surface to pick from
    # boost surf by **3 to enhance probability of picks at peaks (value chosen by testing on synthetics)
    surf = (m.lam1/m.lam2)**3
    dlag = mlags[1] - mlags[0]

    # apply polar density correction
    # not sure about this
    # density = rho(m.lags,dlag)
    # surf = surf / density

    # normalise (probs must add to 1)
    surf = surf / surf.sum()

    # pick fast and tlag from surf
    probs = surf.ravel()
    picks = np.random.choice(probs.size,size=kwargs['n'],replace=True,p=probs)
    idx = np.unravel_index(picks,surf.shape)

    # generate bootstrap sample measurements    
    bslist = [ EigenM(bs,lags=mlags,degs=mdegs) for bs in \
                [ bs_pair(pair,fast,lag) for fast,lag \
                    in zip(m.degs[idx],m.lags[idx]) ] ]
    return bslist

def bs_pair(pair,fast,lag,**kwargs):
    """
    Return data with new noise sequence
    """    
    # copy original data
    bs = pair.copy()   
    origang = bs.cmpangs()[0]
    # replace noise sequence
    bs.unsplit(fast,lag)
    bs.rotateto(bs.pol)
    bs.y = core.resample_noise(bs.y)
    bs.rotateto(origang)
    bs.split(fast,lag)
    return bs


# def rho(n,step):
#     """
#     Polar density of measurements
#     """
#     if n == 0:
#         return 1 / (np.pi/4 * step**2)
#     elif n > 0:
#         return 1 / (2 * np.pi * n * step)
#     else:
#         raise Exception('n not valid')
#
# rho = np.vectorize(rho)
    

    
# def boot_std(listM):
#     avg = np.average(np.stack(listM))
#     diffs = [ M - avg for M in listM ]
#     np.transpose(diffs)
#     ( 1 / ( len(listM) - 1 ) ) * [ np.transpose(m - avg)]
    
