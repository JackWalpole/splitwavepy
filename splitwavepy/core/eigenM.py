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


class EigenM(Measure):
    
    """
    Silver and Chan (1991) eigenvalue method measurement.
    
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
        Populates an EigenM instance.
        """        
        #
        # # process input
        # if len(args) == 1 and isinstance(args[0],Pair):
        #     self.data = args[0]
        # else:
        #     self.data = Pair(*args,**kwargs)
        
        # Derive from Measure
        Measure.__init__(self, data, core.eigvalcov, **kwargs)

        # MAKE MEASUREMENT
        gridout = np.asarray(self.gridsearch(**kwargs))
        self.lam1, self.lam2 = gridout[:,:,1].T, gridout[:,:,0].T
        self.maxloc = core.max_idx(self.lam1/self.lam2)
        
        deggrid, laggrid = self._grid_degs_lags()
        
        #
        # # get some measurement attributes
        # # Using signal to noise ratio in 2-D inspired by 3-D treatment of:
        # # Jackson, Mason, and Greenhalgh, Geophysics (1991)
        # self.snrsurf = (self.lam1-self.lam2) / (2*self.lam2)
        # maxloc = core.max_idx(self.snrsurf)
        self.fast = deggrid[self.maxloc]
        self.lag  = laggrid[self.maxloc]
        # self.snr = self.snrsurf[maxloc]
        # # get errors    
        if bootstrap is True:
            self.conf95level = self.bootstrap_conf95()
            self.errsurf = self.lam1 / self.lam2
            self.dfast, self.dlag = self.get_errors(surftype='max')
        else:
            self.conf95level = self.F_conf95()
            self.errsurf = self.lam2
            self.dfast, self.dlag = self.get_errors(surftype='min')
        
    
    def vals(self):
        """returns standard test values for this method i.e. lam1/lam2."""
        return self.lam1 / self.lam2  

    def F_conf95(self):
        """Value of lam2 at 95% confidence contour."""
        return core.ftest(self.lam2, self.ndf(), alpha=0.05)
        
    def bootstrap_conf95(self, **kwargs):
        """Return lam2 value at 95% confidence level"""
        lam1, lam2 = np.asarray(self._bootstrap_loop(**kwargs))
        return np.percentile(lam1/lam2, 2.5)
        
    # Sandvol and Hearn Bootstrapping
        
    def bootstrap_sandh(self, **kwargs):
        """
        Bootstrap
        """
        
        if 'n' not in kwargs: raise Exception('number of bootstrap iterations *n* required, e.g., n=50')
        
        deggrid, laggrid = self._grid_degs_lags()
        
        resultslist = []
        
        # generate bootstrap sample measurements
        for bs in ( self._renoise(**kwargs) for x in range(kwargs['n']) ):
            gridout = np.asarray(bs.gridsearch(**kwargs))
            bs.lam1, bs.lam2 = gridout[:,:,1].T, gridout[:,:,0].T
            bs.maxloc = core.max_idx(bs.lam1/bs.lam2)
            bs.fast = deggrid[bs.maxloc]
            bs.lag  = laggrid[bs.maxloc]
            resultslist.append((bs.fast, bs.lag))
            
        return resultslist
    
    # def bootstrap_sandh
        
    # auto null classification  
    
    # def ni(self):
    #     """
    #     development.
    #     measure of self-similarity in measurements at 90 degree shift in fast direction
    #     """
    #     fastprof = self.fastprofile()
    #     halfway = int(self.degs.shape[1]/2)
    #     diff = fastprof - np.roll(fastprof,halfway)
    #     mult = fastprof * np.roll(fastprof,halfway)
    #     sumdiffsq = np.sum(diff**2)
    #     summult = np.sum(mult)
    #     return summult/sumdiffsq     
    
    # def bootstrap(self, **kwargs):
    #     if 'n' not in kwargs: kwargs['n'] = 50
    #     bslist = [ EigenM(bs, **self.kwargs) for bs in \
    #                 [ self._bootstrap_sample() for ii in range(kwargs['n']) ] ]
    #     return bslist
        
    # def conf_95(self, **kwargs):
    #     """Value of lam2 at 95% confidence contour."""
    #     if 'n' not in kwargs: kwargs['n'] = 500
    #     bslist = [ bs.unsplit(self.fast, self.lag).eigvalcov() for bs in \
    #                 [ self._bootstrap_sample() for ii in range(kwargs['n']) ] ]
    #     lam2s = np.asarray(bslist)[:,0]
    #     # return lam2s
    #     return core.val_at_alpha(lam2s, 0.975)
    
    # Plotting
    
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
           kwargs['vals'] = self.lam1 / self.lam2
           kwargs['title'] = r'$\lambda_1 / \lambda_2$'
        
        self._plot(**kwargs)
    


        
    # def plot_profiles(self,**kwargs):
    #     # Error analysis
    #     fig,ax = plt.subplots(2)
    #     ax0 = plt.subplot(121)
    #     ax1 = plt.subplot(122)
    #
    #     ax0.plot(self.degs[0,:],self.fastprofile())
    #     ax0.axvline(self.fast)
    #     ax0.axvline(self.fast-2*self.dfast,alpha=0.5)
    #     ax0.axvline(self.fast+2*self.dfast,alpha=0.5)
    #     ax0.set_title('fast direction')
    #
    #     ax1.plot(self.lags[:,0],self.lagprofile())
    #     ax1.axvline(self.lag)
    #     ax1.axvline(self.lag-2*self.dlag,alpha=0.5)
    #     ax1.axvline(self.lag+2*self.dlag,alpha=0.5)
    #     ax1.set_title('lag direction')
    #
    #     plt.show()
        





        

        
