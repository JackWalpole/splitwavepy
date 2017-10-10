# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .eigenM import EigenM

import numpy as np
# import matplotlib.pyplot as plt

class Stack:

    def __init__(self,*args,**kwargs):
        """
        A collection (or stack) of EigenM measurements.
        
        To initate a Stack of measurements stored in a list
        
        >>> S = sw.Stack(listM, weights=weights)
    
        To perform lam1 / lam2 stacking
        
        >>> S.stack()
        
        (Jack Walpole: "This is my own personal preference 
                        as lam1 / lam2 surfaces are self-normalised 
                        and better measurements are naturally
                        more heavily weighted.")
    
        To stack lam2 surfaces following Wolfe and Silver (1998)
        
        >>> S.wolfe_silver()
    
        Optionally set keywords to extend the Wofe and Silver method following Restivo and Helffrich (1999)
        -- snr = True -- signal to noise ratio weights added
        -- azi = True -- azimuthal bin weighting added (geometry information required)
        
        >>> S.wolfe_silver(snr=True, azi=True)

        ----------------------------------
        Positional Arguments
        ----------------------------------

        list of EigenM instances
        
        ----------------------------------
        Keyword Arguments
        ----------------------------------

        weights = numpy array of weights (order and size must match list)

        """

        if not isinstance(args[0],list):
            raise Exception('requires a list of EigenM measurements as the first positional argument')
            
        # 
        self.listM = args[0]
        
        # get degs and tlags from first M
        self.degs = self.listM[0].degs
        self.tlags = self.listM[0].tlags
        
        if 'weights' in kwargs:
            self.weights = kwargs['weights']

    def wolfe_silver(self,**kwargs):
        """
        Return stack using method of Wolfe and Silver (1998).
        
        This method stacks the lambda2 surfaces, 
        pre-normalises surfaces so that minimum lambda2 = 1.
        
        Optionally set keywords to extend the Wofe and Silver method following Restivo and Helffrich (1999)
        -- snr = True -- signal to noise ratio weights added
        -- baz = True -- back-azimuthal bin weighting added (geometry information required)
        >>> S.wolfe_silver(snr=True, azi=True)
        """
        
        listS = [ M.lam2 / np.min(M.lam2) for M in self.listM ]
        
        if kwargs['snr'] is True:
            # weight by signal to noise ratio
            # need sigmoid function with min to max ranging from 1 to 21.
            raise Exception('not yet supported')
                  
        if kwargs['baz'] is True:
            # weight by backazimuthal density coverage
            raise Exception('not yet supported')
            
        # keep track of degrees of freedom    
        # ndf = np.sum([M.ndf for M in self.listM])
        
        return _stack(listS,**kwargs)
        

    def stack(self,**kwargs):
        """
        Return stack of lam1 / lam2.
        """
        
        listS = [ M.lam1 / M.lam2 for M in self.listM ]
        
        return _stack(listS,**kwargs)
   

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
    return np.average(stack,axis=0,**kwargs)


