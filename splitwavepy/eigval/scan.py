# -*- coding: utf-8 -*-
"""
Search window parameter space
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core.pair import Pair
from ..core.window import Window
from . import EigenM

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# class Search:
#     """
#     Class to trial different analysis windows
#     """
#
#     def __init__(self,pair,widths,offsets,*args,**kwargs):
#
#         self.pair = pair
#         self.widths = widths
#         self.offsets = offsets
        
    
def window_scan(pair,widths,offsets,**kwargs):
    """
    Returns a list of EigenM measurements.
    
    trials all combinations of *offsets* and *widths* windows
    """    
    
    if not isinstance(pair,Pair):
        raise Exception('pair must be a Pair')
        
    if not isinstance(offsets,np.ndarray):
        raise Exception('offsets must be a numpy.array')
    
    if not isinstance(widths):
        raise Exception('widths must be a numpy.ndarray')
    
    # populate a list of measurements
    listM=[]
    for w in widths:
        for o in offsets:
            listM.append(sw.EigenM(pair,window=sw.Window(w,o,**kwargs)))
    
    return listM
        

# def frequency_scan():
    

# def receiver_correction_scan():
    

# def source_correction_scan():  


def plot(listM,x,y):
    """
    Plot error surfaces in list
    """
    fig = plt.figure(figsize=(12,6)) 
    nwid = widths.size
    noff = offsets.size
    gs = gridspec.GridSpec(nwid,noff)
    v = np.linspace(0, 50, 26, endpoint=True)

    for ii in range(nwid):
        for jj in range(noff):
            ax = plt.subplot(gs[ii,jj])
            m = listM[ii*noff+jj]
            ax.contourf(m.tlags,m.degs,m.lam1/m.lam2,v,cmap='magma',extend='max')

    plt.show()
    
