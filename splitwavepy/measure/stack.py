# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# from .eigenM import EigenM

import numpy as np
import matplotlib.pyplot as plt

class Stack:

    def __init__(self,*args,**kwargs):
        """
        A collection (or stack) of EigenM measurements.
        
        args
        -----
        
        listM      list, a list of EigenM objects (all on same grid)
        
        keyword arguments
        -----------------
        
        weights     numpy array, user-defined weights (order and size must match list)
               
        example
        --------
      
        >>> S = sw.Stack(listM)           # initiate a Stack of EigenM objects
        >>> l1l2stk = S.stack()           # Lam1 / Lam2 stack
        >>> WSstk = S.wolfe_silver()      # Stack using Wolfe and Silver method
        >>> RHstk = S.restivo_helffrich() # Stack using Restivo and Helffrich method
        """

        if not isinstance(args[0], list):
            raise Exception('requires a list of EigenM measurements as the first positional argument')
            
        # get info
        self.listM = args[0]
        self.degs = self.listM[0].degs
        self.lags = self.listM[0].lags
        
        # check all have the same grids
        for M in self.listM:
            if np.any(M.degs != self.degs):
                raise Exception('Inconsistent degs grid found, all surface must have same grids.')
            if np.any(M.lags != self.lags):
                raise Exception('Inconsistent lags grid found, all surfaces must have same grids.')
        
        # if weights provided check it is the right size
        if 'weights' in kwargs:
            if not isinstance(kwargs['weights'], np.ndarray):
                raise TypeError('weights should be a numpy array')
            if len(self.listM) != kwargs['weights'].size:
                raise Exception('weights array size must equal listM length')
            self.weights = kwargs['weights'] # save weights for easy recall

    def wolfe_silver(self, **kwargs):
        """
        Return stack using method of Wolfe and Silver (1998).

        This method stacks the lambda2 surfaces,
        pre-normalises surfaces so that minimum lambda2 = 1.
        """

        listS = [ M.lam2 / np.min(M.lam2) for M in self.listM ]
        return _stack(listS, **kwargs)

    def restivo_helffrich(self, **kwargs):
        """
        Return stack using method of Restivo and Helffrich (1999).
        
        This method is similar to the Wolfe and Silver method except 
        error surfaces are weighted by their signal to noise ratio.        
        """
        
        listS = [ M.lam2 / np.min(M.lam2) for M in self.listM ]

        # weight by signal to noise ratio
        weights = np.asarray([ M.data.snrRH() for M in self.listM ])
     
        # should apply sigmoid (?) function to weights with min to max ranging from 1 to 21 to be consistent with original paper.  Note: sheba does not bother with this.

        # if 'baz' in kwargs and kwargs['baz'] is True:
        #     # weight by backazimuthal density coverage
        #     raise Exception('not yet supported')

        # if user specified weights then apply these on top:
        if 'weights' in kwargs:
            weights = weights * kwargs['weights']
        
        return _stack(listS, weights=weights)

    def stack(self,**kwargs):
        """
        Return stack of lam1 / lam2 with optional user-defined weights.
        """        
        listS = [ (M.lam1-M.lam2) / M.lam2 for M in self.listM ]        
        return _stack(listS, **kwargs)
   
    def stackpdf(self,**kwargs):
        """
        Return stack of lam1 / lam2 with optional user-defined weights.
        """        
        listS = [ ( (M.lam1-M.lam2) / M.lam2) / np.sum( (M.lam1-M.lam2) / M.lam2) for M in self.listM ]        
        return _stack(listS, **kwargs)

# basic stacking routine
def _stack(listSurfaces,**kwargs):
    """
    Average listSurfaces, optionally using weights.

    **kargs

    weights = numpy array of weights
    ----------------------------------
    """

    # put everything into one giant numpy array
    stack = np.stack(listSurfaces)

    # perform the averaging
    # note np.average takes an optional weights parameter so automatically handles weigths.
    return np.average(stack, axis=0, **kwargs)


