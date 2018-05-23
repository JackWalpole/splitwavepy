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

class XC(Measure):
    
    """
    e.g. Ando and Bowman (1987) rotation correlation method.
    
    Calculates Pearson r.
    
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
    """
    
    def __init__(self, data, bootstrap=False, **kwargs):
               
        # process input
        if 'pol' in kwargs: 
            # print('Warning pol not used.')
            del kwargs['pol']
        
        # Derive from Measure
        # Measure.__init__(self, data, core.crosscorr, **kwargs)
        Measure.__init__(self, data, core.pearson, **kwargs)

        # MAKE MEASUREMENT
        gridout = np.asarray(self.gridsearch(self.func, **kwargs))
        self.xc = np.abs(gridout[:,:,0].T)
        maxloc = core.max_idx(self.xc)

        #
        deggrid, laggrid = self._grid_degs_lags()
        self.fast = deggrid[maxloc]
        self.lag  = laggrid[maxloc]

        # # get errors
        if bootstrap is True:
            self.pdf = self.estimate_pdf(**kwargs)
            self.conf95level = self._pdf_conf95(self.pdf)
            self.errsurf = self.pdf
            self.dfast, self.dlag = self.get_errors(surftype='max')
        else:
            self.conf95level = self.F_conf95()
            self.errsurf = self.xc
            self.dfast, self.dlag = self.get_errors(surftype='max')

    def vals(self):
        return self.xc
    
    def F_conf95(self):
        """
        Uses a Fisher transform to calculate confidence interval
        returns (positive) xc value at 95% confidence interval
        Standard Error, SE = 1 / ( sqrt(n-3))
        """    
        # max xc value
        xcmax = np.max(self.xc)
        # Fisher Transform
        f = math.atanh(xcmax)
        std = 1 / math.sqrt((self.ndf()-3))
        z95 = 1.96
        # reverse transform
        xc95 = math.tanh(f-z95*std) 
        return xc95
        
    # def bootstrap_conf95(self, **kwargs):
    #     """Return lam2 value at 95% confidence level"""
    #     xc = self._bootstrap_loop(**kwargs)
    #     return np.percentile(xc, 2.5)
        
    def _bootstrap_prep(self):
        x, y = self.fastdata_corr().chopdata()
        return x, y   
        
    def _bootstrap_stat(self, x, y):
        xc = self.func(x, y)
        return math.fabs(xc)
        
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
    



        
