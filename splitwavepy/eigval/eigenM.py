# -*- coding: utf-8 -*-
"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
from ..core.pair import Pair
from ..core.window import Window
from . import eigval
from ..plotting import plot
from ..core import io

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class EigenM:
    
    """
    Silver and Chan (1991) eigenvalue method measurement.
    
    args:
    None = create synthetic
    Pair = Measure splitting on Pair object
    x, y = Measure splitting on traces x, and y.
    
    kwargs:
    
    window = Window
    
    lags -- None | defaults: minlag=0, maxlag=win.width/10, nlags=30. 
         -- tuple = (maxlag,)  
         -- tuple = (maxlag,Nlags) 
         -- tuple = (minlag,maxlag,Nlags)
    
    degs -- None | defaults: Ndegs=60
         -- int = Ndegs
    
    rcvcorr = (fast,tlag) | tuple | Receiver Correction
    srccorr = (fast,tlag) | tuple | Source Correction
    
    kwargs for synthetic generation:
    fast = 0.    | float
    tlag = 0.    | float
    pol = 0.     | float
    noise = 0.03 | float
        
    """
    
    def __init__(self,*args,**kwargs):
        """
        Populates an EigenM instance.
        """
        
        # process input
        if len(args) == 1 and isinstance(args[0],Pair):
            self.data = args[0]
        else:
            self.data = Pair(*args,**kwargs)
        
        # convert times to nsamples
        self.delta = self.data.delta 
        
        # window
        default_width = self.data.x.size / 2
        default_offset = 0
        default_tukey = None
        
        if 'window' not in kwargs:
            self.window = Window(default_width,default_offset,default_tukey)
        else:
            if not isinstance(kwargs['window'],Window):
                raise TypeError('window must be a Window type')
            else:
                self.window = kwargs['window']
        
        # lags        
        default_minlag = 0
        default_maxlag = self.window.width / 10
        default_nlags  = 30
        
        if 'lags' not in kwargs:
            lags = np.linspace(default_minlag,default_maxlag,default_nlags)
        else:
            if not isinstance(kwargs['lags'],tuple):
                raise TypeError('lags must be a tuple')
            elif len(kwargs['lags']) == 1:
                lags = np.linspace(default_minlag,kwargs['lags'],default_nlags)
            elif len(kwargs['lags']) == 2:
                lags = np.linspace(default_minlag,*kwargs['lags'])
            elif len(kwargs['lags']) == 3:
                lags = np.linspace(*kwargs['lags'])
            else:
                raise Exception('Too much info in lags keyword')
                
            # round to nearest 2
            samplelags = 2 * np.rint( 0.5 * lags / self.delta )
            # a bit messy to reuse kwargs['lags'] in this way but eigval expects this name 
            kwargs['lags'] = np.unique(samplelags).astype(int)
            
        # degs
        default_mindeg = 0
        default_maxdeg = 180
        default_ndegs = 60
        
        if 'degs' not in kwargs:
            degs = np.linspace(default_mindeg,default_maxdeg,default_ndegs,endpoint=False)
        else:
            if not isinstance(kwargs['degs'],int):
                raise TypeError('degs must be an integer')
            degs = np.linspace(default_mindeg,default_maxdeg,kwargs['degs'],endpoint=False)
                
        # receiver correction            
        if ('rcvcorr' in kwargs):
            if not isinstance(kwargs['rcvcorr'],tuple): raise TypeError('rcvcorr must be tuple')
            if len(kwargs['rcvcorr']) != 2: raise Exception('rcvcorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['rcvcorr']
            samps = core.time2samps(lag,self.delta,'even')
            kwargs['rcvcorr'] = (deg, samps)
            self.rcvcorr = ( deg, samps * self.delta)
        
        # source correction                  
        if ('srccorr' in kwargs):
            if not isinstance(kwargs['srccorr'],tuple): raise TypeError('srccorr must be tuple')
            if len(kwargs['srccorr']) != 2: raise Exception('srccorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['srccorr']
            samps = core.time2samps(lag,self.delta,'even')
            kwargs['srccorr'] = (deg, samps)
            self.srccorr = (deg, samps * self.delta)

            
        # ensure trace1 at zero angle
        self.data.rotateto(0)        
        
        # grid search splitting
        self.degs, self.samplags, self.lam1, self.lam2, self.window = eigval.grideigval(self.data.x,self.data.y,**kwargs)
        # convert sample lags to meaningful time lags
        self.lags = self.samplags * self.delta
                
        # get some measurement attributes
        # uses ratio lam1/lam2 to find optimal fast and lag parameters
        maxloc = core.max_idx((self.lam1-self.lam2)/self.lam2)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        
        # correct the data     
        self.data_corr = self.data.copy()
        # rcv side      
        if 'rcvcorr' in kwargs:
            self.data_corr.unsplit(*kwargs['rcvcorr'])    
        # target layer
        self.data_corr.unsplit(self.fast,self.lag)  
        # src side
        if 'srccorr' in kwargs:
            self.data_corr.unsplit(*kwargs['srccorr'])
        
        # recover source polarisation
        self.srcpol = self.data_corr.pca()
        
        # estimate signal to noise   
        self.snr = np.max((self.lam1-self.lam2)/(self.lam2))
        
        # error estimations
        self.fdfast, self.fdlag = self.f_errors()
        
        # other stuff
        self.units = self.data.units
    
    def srcpoldata(self):
        return Pair(*core.rotate(self.data.x,self.data.y,self.srcpol))
        
    def srcpoldata_corr(self):
        # data_corr = self.data_corr
        return Pair(*core.rotate(self.data_corr.x,self.data_corr.y,self.srcpol))
        
    def fastslowdata(self):
        return Pair(*core.rotate(self.data.x,self.data.y,self.fast))
        
    def fastslowdata_corr(self):
        return Pair(*core.rotate(self.data_corr.x,self.data_corr.y,self.fast))
    
    def snrRH(self):
        """Restivo and Helffrich (1999) signal to noise ratio"""
        d = self.srcpoldata_corr()
        return core.snrRH(*core.chop(d.x,d.y,window=self.window))
        
    # F-test utilities
    
    def ndf(self):
        """Number of degrees of freedom."""
        d = self.srcpoldata_corr()
        return eigval.ndf(d.y,window=self.window)
    
    def lam2_95(self):
        """Value of lam2 at 95% confidence contour."""
        return eigval.ftest(self.lam2,self.ndf(),alpha=0.05)
        
    def f_errors(self):
        """
        Return dfast and dtlag.

        These errors correspond to one sigma in the parameter estimate.

        Calculated by taking a quarter of the width of 95% confidence region (found using F-test).
        """

        # search interval steps
        lag_step = self.lags[1,0]-self.lags[0,0]
        fast_step = self.degs[0,1]-self.degs[0,0]

        # Find nodes where we fall within the 95% confidence region
        confbool = self.lam2 <= self.lam2_95()

        # tlag error
        lagbool = confbool.any(axis=1)
        # last true value - first true value
        truth = np.where(lagbool)[0]
        fdlag = (truth[-1] - truth[0] + 1) * lag_step * 0.25

        # fast error
        fastbool = confbool.any(axis=0)
        # trickier to handle due to cyclicity of angles
        # search for the longest continuous line of False values
        cyclic = np.hstack((fastbool,fastbool))
        lengthFalse = np.diff(np.where(cyclic)).max()
        # shortest line that contains ALL true values is then:
        lengthTrue = fastbool.size - lengthFalse
        fdfast = lengthTrue * fast_step * 0.25

        # return
        return fdfast, fdlag
        
        
    # "squashed" profiles
    
    def fastprofile(self):
        return np.sum((self.lam1-self.lam2)/self.lam2, axis=0)
        
    def lagprofile(self):
        return np.sum((self.lam1-self.la2)/self.lam2, axis=1)
    
    # auto null classification  
    
    def ni(self):
        """
        development.
        measure of self-similarity in measurements at 90 degree shift in fast direction
        """
        fastprof = self.fastprofile()
        halfway = int(self.degs.shape[1]/2)
        diff = fastprof - np.roll(fastprof,halfway)
        mult = fastprof * np.roll(fastprof,halfway)
        sumdiffsq = np.sum(diff**2)
        summult = np.sum(mult)
        return summult/sumdiffsq
    
    # Output
    
    def report(self):
        """
        Print out a summary of the result.
        """
        
        # pformat = 

    def save(self,filename):
        """
        Save Measurement for future referral
        """
        io.save(self,filename)
        
    
    # Plotting

    def plot(self,**kwargs):
        
        if 'mode' not in kwargs:
            kwargs['mode'] = 'all'
        
        # plot only surface    
        if kwargs['mode'] == 'surf':

            if 'vals' not in kwargs:
                kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2        
            plot.surf(self,**kwargs)

        # else plot everything   
        else:    
            fig = plt.figure(figsize=(12,6)) 
            gs = gridspec.GridSpec(2, 3,
                               width_ratios=[1,1,2]
                               )
    
            ax0 = plt.subplot(gs[0,0])
            ax1 = plt.subplot(gs[0,1])
            ax2 = plt.subplot(gs[1,0])
            ax3 = plt.subplot(gs[1,1])
            ax4 = plt.subplot(gs[:,2])
            
            # d1 = self.data.copy()
            # d1.chop(self.window)
            # d2 = self.data_corr.copy()
            # d2.chop(self.window)
            d1 = self.srcpoldata()
            d1.chop(self.window)
            d2 = self.srcpoldata_corr()
            d2.chop(self.window)
            
            # get axis scaling
            lim = np.abs(np.hstack((d1.xy(),d2.xy()))).max() * 1.1
            ylim = [-lim,lim]
    
            # original
            plot.trace(d1.x,d1.y,time=d1.t(),ax=ax0,ylim=ylim)
            plot.particle(d1.x,d1.y,ax=ax1,lim=ylim)
    
            # corrected
            plot.trace(d2.x,d2.y,time=d2.t(),ax=ax2,ylim=ylim)
            plot.particle(d2.x,d2.y,ax=ax3,lim=ylim)
    
            # error surface
            if 'vals' not in kwargs:
                kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2
                kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
        
            plot.surf(self,ax=ax4,**kwargs)
        
            # neaten
            plt.tight_layout()

        plt.show()



    # Comparison
    
    def __eq__(self, other) :
        # check same class
        if self.__class__ != other.__class__: return False
        # check same keys
        if self.__dict__.keys() != other.__dict__.keys(): return False
        # check same values
        for key in self.__dict__.keys():
            if np.all( self.__dict__[key] != other.__dict__[key]): return False
        # if reached here then the same
        return True
        



        
