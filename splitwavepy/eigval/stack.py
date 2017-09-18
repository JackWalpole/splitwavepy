"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# from ..core import core
# from ..core import pair
# from ..core import window
from .eigenM import EigenM

import numpy as np
import matplotlib.pyplot as plt

class Stack:
    
    def __init__(self,*args,**kwargs):
        """
        Stack list of EigenM measurements.
        
        ----------------------------------
        *args 
        
        list of EigenM instances
        ----------------------------------
        **kargs
        
        weights = numpy array of weights (order must match list)
        ----------------------------------
        """
        
        if not isinstance(args[0],list):
            raise Exception('requires a list of EigenM measurements as the first positional argument')
        
        # list of measurements  
        self.list = args[0]
        
        # numpy array of weights
        if (len(args) > 1 and 
                isinstance(args[1],np.ndarray) and
                    args[1].size == len(self.list)):            
            self.weights = args[1]
        else:
            self.weights = np.ones(len(self.list))
        
        # stacking routine
        def stack(self,**kwargs):
            """
            Return numpy array of stacked lam1/lam2 values
            
            **kwargs
            """  
            stack = np.zeros(listM[0].degs.shape)            
            for ii in range(len(self.list)):
                m = self.list[ii]         
                stack += self.weight[ii] * (m.lam1/m.lam2)
            # 
            return stack / np.sum(self.weights)


