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
# from ..plotting import plot
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
    fast = 0.      | float
    tlag = 0.      | float
    pol = 0.       | float
    noise = 0.001  | float
        
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
        self.units = self.data.units
        
        # window default parameters
        width = self.data.x.size / 2
        offset = 0
        tukey = None
        
        if 'window' not in kwargs:
            self.window = Window( width, offset, tukey)
        else:
            if not isinstance(kwargs['window'],Window):
                raise TypeError('window must be a Window type')
            else:
                self.window = kwargs['window']
        
        # lags        
        minlag = 0
        maxlag = self.window.width / 10
        nlags  = 30
        
        if 'lags' not in kwargs:
            lags = np.linspace( minlag, maxlag, nlags)
        else:
            if not isinstance(kwargs['lags'],tuple):
                raise TypeError('lags must be a tuple')
            elif len(kwargs['lags']) == 1:
                lags = np.linspace( minlag, kwargs['lags'], nlags)
            elif len(kwargs['lags']) == 2:
                lags = np.linspace( minlag, *kwargs['lags'])
            elif len(kwargs['lags']) == 3:
                lags = np.linspace( *kwargs['lags'])
            else:
                raise Exception('Too much info in lags keyword')
            
            # lags must be even
            # a bit messy to reuse kwargs['lags'] in this way but eigval expects this name 
            kwargs['lags'] = np.unique( core.time2samps( lags, delta, mode='even')).astype(int)
            
        # degs
        mindeg = 0
        maxdeg = 180
        ndegs = 60
        
        if 'degs' not in kwargs:
            degs = np.linspace( mindeg, maxdeg, ndegs, endpoint=False)
        else:
            if not isinstance(kwargs['degs'],int):
                raise TypeError('degs must be an integer')
            degs = np.linspace( mindeg, maxdeg, kwargs['degs'], endpoint=False)
                
        # receiver correction            
        if ('rcvcorr' in kwargs):
            if not isinstance(kwargs['rcvcorr'],tuple): raise TypeError('rcvcorr must be tuple')
            if len(kwargs['rcvcorr']) != 2: raise Exception('rcvcorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['rcvcorr']
            samps = core.time2samps( lag,self.delta, 'even')
            kwargs['rcvcorr'] = ( deg, samps)
            self.rcvcorr = ( deg, samps * self.delta)
        
        # source correction                  
        if ('srccorr' in kwargs):
            if not isinstance(kwargs['srccorr'],tuple): raise TypeError('srccorr must be tuple')
            if len(kwargs['srccorr']) != 2: raise Exception('srccorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['srccorr']
            samps = core.time2samps( lag, self.delta, 'even')
            kwargs['srccorr'] = ( deg, samps)
            self.srccorr = ( deg, samps * self.delta)

            
        # ensure trace1 at zero angle
        self.data.rotateto(0)        
        
        # grid search splitting
        self.degs, self.samplags, self.lam1, self.lam2, self.window = eigval.grideigval(self.data.x,self.data.y,**kwargs)
        # convert sample lags to meaningful time lags
        self.lags = self.samplags * self.delta
                
        # get some measurement attributes
        # maximise lam1/lam2
        maxloc = core.max_idx(self.lam1/self.lam2)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        # estimate signal to noise level (lam1-lam2)/lam2
        self.snr = (self.lam1[maxloc]-self.lam2[maxloc])/(self.lam2[maxloc])
        
        # data corrections
        self.rcvcorr = None
        if 'rcvcorr' in kwargs: self.rcvcorr = kwargs['rcvcorr']
        self.srccorr = None
        if 'srccorr' in kwargs: self.rcvcorr = kwargs['srccorr'] 
        

    # methods
    
    def srcpol(self):
        # recover source polarisation
        return self.data_corr().pol()
        
    def snrRH(self):
        """Restivo and Helffrich (1999) signal to noise ratio"""
        d = self.srcpoldata_corr()
        x,y = d.chop()
        return core.snrRH(x,y)
                
    # possibly useful rotations
    
    def data_corr(self):        
        # copy data     
        data_corr = self.data.copy()
        # rcv side correction     
        if self.rcvcorr is not None:
            data_corr.unsplit(*self.rcvcorr)    
        # target layer correction
        data_corr.unsplit(self.fast,self.lag)  
        # src side correction
        if self.srccorr is not None:
            data_corr.unsplit(*self.srccorr)
        return data_corr

    def srcpoldata(self):
        srcpoldata = self.data.copy()
        srcpoldata.rotateto(self.srcpol())
        return srcpoldata
        
    def srcpoldata_corr(self):
        srcpoldata_corr = self.data_corr()        
        srcpoldata_corr.rotateto(self.srcpol())
        return srcpoldata_corr
        
    def fastdata(self):
        fastdata = self.data.copy()
        fastdata.rotateto(self.fast)
        return fastdata

    def fastdata_corr(self):
        fastdata_corr = self.data_corr()
        fastdata_corr.rotateto(self.fast)
        return fastdata_corr
            
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
          
        # parse kwargs      
        
        if 'vals' not in kwargs:
            kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2        
 
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
        
        # sata to plot
        d1 = self.data.chop()
        # d1 = self.fastdata().chop()
        d2 = self.data_corr().chop()
        # d3 = self.srcpoldata_corr().chop()
        
        # get axis scaling
        lim = np.abs(np.hstack((d1.data(),d2.data()))).max() * 1.1
        ylim = [-lim,lim]

        # original
        d1._ptr(ax0,ylim=ylim,**kwargs)
        d1._ppm(ax1,lims=ylim,**kwargs)
        # corrected
        d2._ptr(ax2,ylim=ylim,**kwargs)
        d2._ppm(ax3,lims=ylim,**kwargs)

        # error surface
        if 'vals' not in kwargs:
            kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2
            kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'

        self._psurf(ax4,**kwargs)
    
        # neaten
        plt.tight_layout()
        plt.show()
        
    def _psurf(self,ax,**kwargs):
        """
        Plot an error surface.
    
        **kwargs
        - cmap = 'magma'
        - vals = (M.lam1-M.lam2) / M.lam2
        - ax = None (creates new)
        """
    
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'magma'
    
        if 'vals' not in kwargs:
            kwargs['vals'] = (self.lam1-self.lam2) / self.lam2
            
        # error surface
        cax = ax.contourf(self.lags,self.degs,kwargs['vals'],26,cmap=kwargs['cmap'])
        cbar = plt.colorbar(cax)
        ax.set_yticks(np.linspace(-90,90,6,endpoint=False))
        ax.set_ylabel(r'Fast Direction ($^\circ$)')
        ax.set_xlabel('Delay Time (' + self.units + ')')
        
        # # marker
        # ax.errorbar(M.lag,M.fast,xerr=M.fdlag,yerr=M.fdfast,fmt='o')
        #
        # # confidence region
        # ax.contour(M.lags,M.degs,M.lam2,levels=[M.lam2_95()])

        ax.set_xlim([self.lags[0,0], self.lags[-1,0]])
        ax.set_ylim([self.degs[0,0], self.degs[0,-1]])
    
        # optional title
        if 'title' in kwargs:
            ax.set_title(kwargs['title'])

        return ax

#     def plot(self,**kwargs):
#
#
#         if 'cmplabels' not in kwargs:
#             kwargs['cmplabels'] = self.data.cmplabels
#
#         if 'mode' not in kwargs:
#             kwargs['mode'] = 'all'
#
#         # plot only surface
#         if kwargs['mode'] == 'surf':
#
#             if 'vals' not in kwargs:
#                 kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2
#             plot.surf(self,**kwargs)
#
#         # else plot everything
#         else:
#             fig = plt.figure(figsize=(12,6))
#             gs = gridspec.GridSpec(2, 3,
#                                width_ratios=[1,1,2]
#                                )
#
#             ax0 = plt.subplot(gs[0,0])
#             ax1 = plt.subplot(gs[0,1])
#             ax2 = plt.subplot(gs[1,0])
#             ax3 = plt.subplot(gs[1,1])
#             ax4 = plt.subplot(gs[:,2])
#
#             d1 = self.srcpoldata()
#             d2 = self.srcpoldata_corr()
#
#             d1d = d1.chop()
#             d1t = d1.chopt()
#             d2d = d2.chop()
#             d2t = d2.chopt()
#
#             # get axis scaling
#             lim = np.abs(np.hstack((d1d,d2d))).max() * 1.1
#             ylim = [-lim,lim]
#
# ## TODO fix times so they are equal at centre sample
#
#             # original
#             plot.trace(d1d[0],d1d[1],time=d1t,ax=ax0,ylim=ylim,**kwargs)
#             plot.particle(d1d[0],d1d[1],ax=ax1,lim=ylim,**kwargs)
#
#             # corrected
#             plot.trace(d2d[0],d2d[1],time=d2t,ax=ax2,ylim=ylim,**kwargs)
#             plot.particle(d2d[0],d2d[1],ax=ax3,lim=ylim,**kwargs)
#
#             # error surface
#             if 'vals' not in kwargs:
#                 kwargs['vals'] = (self.lam1 - self.lam2) / self.lam2
#                 kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
#
#             plot.surf(self,ax=ax4,**kwargs)
#
#             # neaten
#             plt.tight_layout()
#
#         plt.show()



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
        



        
