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
    
    # grunt work using list comprehension
    listM = [ [ EigenM(pair,window=Window(w,o,**kwargs)) for w in widths ] for o in offsets ]
    
    return listM
        

       
            
    