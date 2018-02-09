# -*- coding: utf-8 -*-
"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
from .measure import Measure

import numpy as np
import math

class XcorrM(Measure):
    
    """
    e.g. Ando and Bowman (1987) rotation correlation method.
    
    Usage:
    m = CrossM(pair)
    
    args:
    None = create synthetic
    Pair = Measure splitting on Pair object
    x, y = Measure splitting on traces x, and y.
    
    kwargs:
    
    name -- string = 'Untitled'
    
    lags -- tuple = (maxlag,)  
         -- tuple = (maxlag,Nlags) 
         -- tuple = (minlag,maxlag,Nlags)
         -- numpy ndarray
    
    degs -- int = degs
         -- numpy ndarray
    
    rcvcorr = (fast,tlag) | tuple | Receiver Correction
    srccorr = (fast,tlag) | tuple | Source Correction
    
    kwargs for synthetic generation:
    fast = 0.      | float
    tlag = 0.      | float
    pol = 0.       | float
    noise = 0.001  | float
        
    """
    
    def __init__(self,*args,**kwargs):
        
        #
        # # process input
        # if len(args) == 1 and isinstance(args[0],Pair):
        #     self.data = args[0]
        # else:
        #     self.data = Pair(*args,**kwargs)
        
        # Derive from Measure
        Measure.__init__(self, *args, **kwargs)

        # MAKE MEASUREMENT
        stuff = np.asarray(self.gridsearch(core.crosscorr,**kwargs))
        self.xc = np.abs(stuff[:,:,0].T)
        maxloc = core.max_idx(self.xc)

        #
        deggrid, laggrid = self._grid_degs_lags()
        self.fast = deggrid[maxloc]
        self.lag  = laggrid[maxloc]

        # # get errors
        self.errsurf = self.xc
        self.dfast, self.dlag = self.get_errors(surftype='max')

        # Name
        self.name = 'Untitled'
        if 'name' in kwargs: self.name = kwargs['name']


    
    def conf_95(self):
        """
        returns (positive) xc value at 95% confidence interval
        """    
        # max xc value
        xcmax = np.max(self.xc)
        # Fisher Transform
        z = math.atanh(xcmax)
        std = math.sqrt( 1 / (self.ndf()-3) )
        z95 = z - 1.96 * std
        xc95 = math.tanh(z95) 
        return xc95
        
            
    
    # Plotting
    
    # Plotting
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
            kwargs['vals'] = self.xc
            kwargs['title'] = 'Correlation Coefficient'
        
        self._plot(**kwargs)
    



        
