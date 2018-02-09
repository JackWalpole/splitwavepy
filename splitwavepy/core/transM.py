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
    
    def __init__(self,*args,**kwargs):
        """
        Populates a TransM instance.
        """
        
        # process input
        if 'pol' not in kwargs: raise Exception('Polarisation must be specified using the pol keyword argument, e.g., pol=30.')
        # self.pol = kwargs['pol']
        
        #
        # # process input
        # if len(args) == 1 and isinstance(args[0],Pair):
        #     self.data = args[0]
        # else:
        #     self.data = Pair(*args,**kwargs)
        
        # Derive from Measure
        Measure.__init__(self, *args, **kwargs)

        # MAKE MEASUREMENT
        stuff = np.asarray(self.gridsearch(core.transenergy, mode='rotpol', **kwargs))
        self.energy1, self.energy2 = stuff[:,:,0].T, stuff[:,:,1].T
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
        self.dfast, self.dlag = self.get_errors(surftype='min')

        # Name
        self.name = 'Untitled'
        if 'name' in kwargs: self.name = kwargs['name']


    def conf_95(self):
        """Value of energy2 at 95% confidence contour."""
        return core.ftest(self.energy2, self.ndf(), alpha=0.05)       
    
    # Plotting
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
            # kwargs['vals'] = (self.energy1 - self.energy2) / self.energy2
            # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
            kwargs['vals'] = self.energy1 / self.energy2
            kwargs['title'] = r'pol energy / trans energy'
        
        self._plot(**kwargs)

        




        
