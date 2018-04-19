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
    
    def __init__(self, data, bootstrap=True, **kwargs):
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
        self.errsurf = self.energy2
        if bootstrap is True:
            self.conf95level = self.bootstrap_conf95()
        else:
            self.conf95level = self.conf_95()
        self.dfast, self.dlag = self.get_errors(surftype='min')

        # Name
        self.name = 'Untitled'
        if 'name' in kwargs: self.name = kwargs['name']


    def conf_95(self):
        """Value of energy2 at 95% confidence contour."""
        return core.ftest(self.energy2, self.ndf(), alpha=0.05)
        
    def bootstrap_conf95(self, **kwargs):
        """Return lam2 value at 95% confidence level"""
        vals = np.asarray(self._bootstrap_loop(**kwargs))
        lam2 = vals[:,1]
        lam1 = vals[:,0]
        return np.percentile(lam2, 97.5)   
    
    # Plotting
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
            # kwargs['vals'] = (self.energy1 - self.energy2) / self.energy2
            # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
            kwargs['vals'] = self.energy1 / self.energy2
            kwargs['title'] = r'pol energy / trans energy'
        
        self._plot(**kwargs)

        




        
