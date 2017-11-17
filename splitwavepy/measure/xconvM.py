# -*- coding: utf-8 -*-
"""
Grid search for params that best model waveform using cross-convolution method.
Menke and Levin (2004)
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core, io
from ..core.pair import Pair
from ..core.window import Window
from .measure import Measure

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os.path

class ConvM(Measure):
    
    """
    Menke and Levin (1991) single layer cross convolution grid search.
    
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
        if 'pol' not in kwargs: raise Exception('Polarisation must be specified, e.g., pol=30.')
        # self.pol = kwargs['pol']
  
        # process input
        if len(args) == 1 and isinstance(args[0],Pair):
            self.data = args[0]
        else:
            self.data = Pair(*args,**kwargs)
        
        # Derive from Measure
        Measure.__init__(self, *args, **kwargs)

        # MAKE MEASUREMENT
        stuff = np.asarray(self.onelayer(core.crossconvmf,**kwargs))
        self.mf = stuff[:,:,0].T
        maxloc = core.min_idx(self.mf)
        
        #
        # # get some measurement attributes
        # # Using signal to noise ratio in 2-D inspired by 3-D treatment of:
        # # Jackson, Mason, and Greenhalgh, Geophysics (1991)
        # self.snrsurf = (self.lam1-self.lam2) / (2*self.lam2)
        # maxloc = core.max_idx(self.snrsurf)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        # self.snr = self.snrsurf[maxloc]
        # # get errors
        self.errsurf = self.lam2
        self.dfast, self.dlag = self.get_errors(surftype='min')
        
        # Name
        self.name = 'Untitled'
        if 'name' in kwargs: self.name = kwargs['name']


    # def conf_95(self):
    #     """Value of lam2 at 95% confidence contour."""
    #     return core.ftest(self.lam2, self.ndf(), alpha=0.05)
    #
    # # auto null classification
    #
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
    
    #
    
    def gridsearchsynth(self, func, **kwargs):
        
        """
        Grid search for splitting parameters applied to data using the function defined in func
        rcvcorr = receiver correction parameters in tuple (fast,lag) 
        srccorr = source correction parameters in tuple (fast,lag) 
        """
        
        # avoid using "dots" in loops for performance
        synth = core.synth        
        # rotate = core.rotate
        # lag = core.lag
        chop = core.chop
        unsplit = core.unsplit
        
        # ensure trace1 at zero angle
        copy = self.data.copy()
        copy.rotateto(0)
        x, y = copy.x, copy.y
        
        # pre-apply receiver correction
        if 'rcvcorr' in kwargs:
            rcvphi, rcvlag = self.__rcvcorr
            x, y = unsplit(x, y, rcvphi, rcvlag)
         
        ######################                  
        # inner loop function
        ######################
    
        # source correction  
        
        if 'srccorr' in kwargs:
            srcphi, srclag = self.__srccorr
            def srccorr(x, y, ang):
                x, y = unsplit(x, y, srcphi-ang, srclag)
                return x, y
        else:
            def srccorr(x, y, ang):
                return x, y
                
        
        # actual inner loop function   
        def getout(x, y, ang, shift):
            # generate synthetics
            mx, my = synth(split=(ang,shift), pol=pol, noise=0, delta=delta)
            mx, my = chop(mx, my, window=self.data.window)
            x, y = srccorr(x, y, ang)
            x, y = chop(x, y, window=self.data.window)
            return func(x, y)

                    
        # Do the grid search
        prerot = [ (rotate(x, y, ang), ang) for ang in self.__degs ]
        
        out = [ [ getout(data[0], data[1], ang, shift) for shift in self.__slags ]
                for (data,ang) in prerot  ]
                               
        return out
    
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
            # kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2
            # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
            kwargs['vals'] = self.lam1 / self.lam2
            kwargs['title'] = r'Misfit'
        
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
        





        

        
