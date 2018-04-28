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
    
    def __init__(self, data, bootstrap=False, **kwargs):
        
        # Derive from Measure
        # Measure.__init__(self, data, core.crosscorr, **kwargs)
        Measure.__init__(self, data, core.pearson, **kwargs)

        # MAKE MEASUREMENT
        gridout = np.asarray(self.gridsearch(**kwargs))
        self.xc = np.abs(gridout[:,:,0].T)
        maxloc = core.max_idx(self.xc)

        #
        deggrid, laggrid = self._grid_degs_lags()
        self.fast = deggrid[maxloc]
        self.lag  = laggrid[maxloc]

        # # get errors
        self.errsurf = self.xc
        if bootstrap is True:
            self.conf95level = self.bootstrap_conf95()
        else:
            self.conf95level = self.conf_95()
        self.dfast, self.dlag = self.get_errors(surftype='max')

    def vals(self):
        return self.xc
    
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
        
    def bootstrap_conf95(self, **kwargs):
        """Return lam2 value at 95% confidence level"""
        xc = np.abs(self._bootstrap_loop(**kwargs))
        return np.percentile(xc, 2.5)
        
    def fisher(self):
        return np.arctanh(self.xc)
    
    # Plotting
    
    # Plotting
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
            # kwargs['vals'] = self.xc
            # kwargs['title'] = 'Correlation Coefficient'
            kwargs['vals'] = self.fisher()
            kwargs['title'] = 'Fisher Transformed Correlation Coefficient'
        
        self._plot(**kwargs)
    



        
