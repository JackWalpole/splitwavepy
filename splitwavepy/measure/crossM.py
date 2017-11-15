# -*- coding: utf-8 -*-
"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core,io
from ..core.pair import Pair
from ..core.window import Window
from .measure import Measure

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# import os.path
import math


class CrossM(Measure):
    
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
        
        # process input
        if len(args) == 1 and isinstance(args[0],Pair):
            self.data = args[0]
        else:
            self.data = Pair(*args,**kwargs)
        
        # Derive from Measure
        Measure.__init__(self, *args, **kwargs)

        # MAKE MEASUREMENT
        stuff = np.asarray(self.gridsearch(core.crosscorr,**kwargs))
        self.xc = np.abs(stuff[:,:,0].T)
        maxloc = core.max_idx(self.xc)

        #
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]

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
    
    def plot(self,**kwargs):
          
        # setup figure and subplots
        fig = plt.figure(figsize=(12,6)) 
        gs = gridspec.GridSpec(2, 3,
                           width_ratios=[1,1,2]
                           )    
        ax0 = plt.subplot(gs[0,0])
        ax1 = plt.subplot(gs[0,1])
        ax2 = plt.subplot(gs[1,0])
        ax3 = plt.subplot(gs[1,1])
        ax4 = plt.subplot(gs[:,2])
        
        # data to plot
        d1 = self.data.chop()
        d1f = self.srcpoldata().chop()
        d2 = self.data_corr().chop()
        d2s = self.srcpoldata_corr().chop()
        
        # flip polarity of slow wave in panel one if opposite to fast
        # d1f.y = d1f.y * np.sign(np.tan(self.srcpol()-self.fast))
        
        # get axis scaling
        lim = np.abs(d2s.data()).max() * 1.1
        ylim = [-lim,lim]

        # original
        d1f._ptr(ax0,ylim=ylim,**kwargs)
        d1._ppm(ax1,lims=ylim,**kwargs)
        # corrected
        d2s._ptr(ax2,ylim=ylim,**kwargs)
        d2._ppm(ax3,lims=ylim,**kwargs)

        # error surface
        if 'vals' not in kwargs:
            kwargs['vals'] = self.xc
            kwargs['title'] = 'Correlation Coefficient'
        
        # add marker and info box by default
        if 'marker' not in kwargs: kwargs['marker'] = True
        if 'info' not in kwargs: kwargs['info'] = True
        if 'conf95' not in kwargs: kwargs['conf95'] = True
        self._psurf(ax4,**kwargs)
        
        # title
        if self.name != 'Untitled':
            plt.suptitle(self.name)
        
        # neaten
        plt.tight_layout()
        plt.show()
        


        
