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
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec


class TransM(Measure):
    
    """
    Silver and Chan (1991) transverse minimisation method.
    
    requires polarisation.
    
    With data:
    
    TransM(data, pol)
    
    For synthetic:
    
    TransM(pol, **kwargs)
    
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
        """
        Populates a TransM instance.
        """
        
        # process input
        if 'pol' not in kwargs: 
            raise Exception('Polarisation must be specified using the pol keyword argument, e.g., pol=30.')

        # Derive from Measure
        Measure.__init__(self, data, core.transenergy, **kwargs)

        # MAKE MEASUREMENT
        gridout = np.asarray(self.gridsearch(mode='rotpol', **kwargs))
        self.energy1, self.energy2 = gridout[:,:,0].T, gridout[:,:,1].T
        maxloc = core.max_idx(self.energy1/self.energy2)
        
        #
        # # get some measurement attributes
        # # Using signal to noise ratio in 2-D inspired by 3-D treatment of:
        # # Jackson, Mason, and Greenhalgh, Geophysics (1991)
        # self.snrsurf = (self.energy1-self.energy2) / (2*self.energy2)
        # maxloc = core.max_idx(self.snrsurf)
        deggrid, laggrid = self._grid_degs_lags()
        self.fast = deggrid[maxloc]
        self.lag  = laggrid[maxloc]
        # self.snr = self.snrsurf[maxloc]
        # get errors
        if bootstrap is True:
            self.conf95level = self.bootstrap_conf95()
            self.errsurf = self.energy1 / self.energy2
            self.dfast, self.dlag = self.get_errors(surftype='max')
        else:
            self.conf95level = self.F_conf95()
            self.errsurf = self.energy2
            self.dfast, self.dlag = self.get_errors(surftype='min')

    def vals(self):
        """returns standard test values for this method i.e. energy1/energy2."""
        return self.energy1 / self.energy2


    def F_conf95(self):
        """Value of energy2 at 95% confidence contour."""
        return core.ftest(self.energy2, self.ndf(), alpha=0.05)
        
    def bootstrap_conf95(self, **kwargs):
        """Return energy1/energy2 value at 95% confidence level"""
        sig1, sig2 = self._bootstrap_loop(**kwargs)
        return np.percentile(sig1/sig2, 2.5)
        
    def _bootstrap_stat(self, x, y):
        sig1, sig2 = self.func(x, y)
        return sig1 / sig2
        
    # def bootstrap_samps(self, **kwargs):
    #     sig1, sig2 = self._bootstrap_loop(**kwargs)
    #     return sig1 / sig2
    
    # Plotting
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
            # kwargs['vals'] = (self.energy1 - self.energy2) / self.energy2
            # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
            kwargs['vals'] = self.energy1 / self.energy2
            kwargs['title'] = r'pol energy / trans energy'
        
        self._plot(**kwargs)

        




        
