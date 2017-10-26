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
from . import eigval

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os.path


class EigenM:
    
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
                
        # LAGS               
        minlag = 0
        maxlag = self.data.wwidth() / 4
        nlags  = 40
        if 'lags' not in kwargs:
            lags = np.linspace( minlag, maxlag, nlags)
        else:
            if isinstance(kwargs['lags'],np.ndarray):
                lags = kwargs['lags']                
            elif isinstance(kwargs['lags'],tuple):                
                if len(kwargs['lags']) == 1:
                    lags = np.linspace( minlag, kwargs['lags'][0], nlags)
                elif len(kwargs['lags']) == 2:
                    lags = np.linspace( minlag,*kwargs['lags'])
                elif len(kwargs['lags']) == 3:
                    lags = np.linspace( *kwargs['lags'])
                else:
                    raise Exception('Can\'t parse lags keyword')                   
            else:
                raise TypeError('lags keyword must be a tuple or numpy array')                
        # convert lags to samples (must be even)
        slags = np.unique( core.time2samps( lags, self.delta, mode='even')).astype(int)
        
        # DEGS
        mindeg = -90
        maxdeg = 90
        degs = 90        
        if 'degs' not in kwargs:
            degs = np.linspace( mindeg, maxdeg, degs, endpoint=False)
        else:
            if isinstance(kwargs['degs'],np.ndarray):
                degs = kwargs['degs']
            elif isinstance(kwargs['degs'],int):
                degs = np.linspace( mindeg, maxdeg, kwargs['degs'], endpoint=False)
            else:
                raise TypeError('degs must be an integer or numpy array') 
        sdegs=degs
        kwargs.pop('degs', None)
                
        # receiver correction 
        self.rcvcorr = None           
        if ('rcvcorr' in kwargs):
            if not isinstance(kwargs['rcvcorr'],tuple): raise TypeError('rcvcorr must be tuple')
            if len(kwargs['rcvcorr']) != 2: raise Exception('rcvcorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['rcvcorr']
            samps = core.time2samps( lag,self.delta, 'even')
            kwargs['rcvcorr'] = ( deg, samps)
            self.rcvcorr = ( deg, samps * self.delta)
        
        # source correction
        self.srccorr = None                  
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

        # MAKE MEASUREMENT
        window = self.data.window
        self.degs, self.samplags, self.lam1, self.lam2 \
        = eigval.grideigval(self.data.x, self.data.y, sdegs, slags, window, **kwargs)
        # convert sample lags to meaningful time lags
        self.lags = self.samplags * self.delta
                
        # get some measurement attributes
        # Using signal to noise ratio in 2-D inspired by 3-D treatment of:
        # Jackson, Mason, and Greenhalgh, Geophysics (1991)
        self.snrsurf = (self.lam1-self.lam2) / (2*self.lam2)
        maxloc = core.max_idx(self.snrsurf)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        self.snr = self.snrsurf[maxloc]
        # get errors
        self.dfast, self.dlag = self.f_errors()
        
        # Name
        self.name = 'Untitled'
        if 'name' in kwargs: self.name = kwargs['name']

    
    # METHODS 
    #---------    
                
    # def report(self,fname=None,**kwargs):
    #     """
    #     Report to stdout or to a file.
    #
    #     keywords
    #     --------
    #
    #     fname    string e.g. 'myfile.txt'. If None will write to stdout.
    #     append   bool   e.g. True. Append to existing file.
    #     header   bool   e.g. False.  Include the header line.
    #     choose   list   e.g. ['fast','lag'].  Choose which attributes (and order) to report.
    #
    #     By default will report to stdout with a header.
    #
    #     If a file name is provided using the keyword *file*
    #     then the code, by default, will write with a header
    #     to a new file, and append without a header to a pre-
    #     existing file.
    #
    #     By default the code will report:
    #     name, fast, lag, dfast, dlag,
    #
    #     choose
    #
    #
    #     """
    #
    #     # by default write to stdout and include the header
    #     header = True
    #     append = False
    #
    #     if fname is not None:
    #         if not isinstance(kwargs['file'],str):
    #             raise TypeError('file name must be a string')
    #         # does file exist?
    #         if os.path.isfile(fname):
    #             # yes -- change defaults
    #             header = False
    #             append = True
    #
    #     # overwrite defaults with keyword arguments
    #     if 'header' in kwargs: header = kwargs['header']
    #     if 'append' in kwargs: append = kwargs['append']
    #
    #     # choose what to report
    #     choose=['name','fast','dfast','lag','dlag','snr','ndf','rcvcorr','srccorr']
    #     # get header line
    #     # get data line
    #
    #     # if file not exist
    #     if append
    #
    #     # if file exists
    #     # exists append

    
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
        srcpoldata.set_labels(['srcpol','residual'])
        return srcpoldata
        
    def srcpoldata_corr(self):
        srcpoldata_corr = self.data_corr()        
        srcpoldata_corr.rotateto(self.srcpol())
        srcpoldata_corr.set_labels(['srcpol','residual'])
        return srcpoldata_corr
        
    def fastdata(self,flipslow=False):
        """Plot fast/slow data."""
        fastdata = self.data.copy()
        fastdata.rotateto(self.fast)
        fastdata.set_labels(['fast','slow'])
        return fastdata

    def fastdata_corr(self,flipslow=False):
        fastdata_corr = self.data_corr()
        fastdata_corr.rotateto(self.fast)
        fastdata_corr.set_labels(['fast','slow'])
        return fastdata_corr
            
    # F-test utilities
    
    def ndf(self):
        """Number of degrees of freedom."""
        d = self.srcpoldata_corr().chop()
        return eigval.ndf(d.y)
    
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
        lengthFalse = np.diff(np.where(cyclic)).max() - 1
        # shortest line that contains ALL true values is then:
        lengthTrue = fastbool.size - lengthFalse
        fdfast = lengthTrue * fast_step * 0.25

        # return
        return fdfast, fdlag
        
        
    # "squashed" profiles
    
    def fastprofile(self):
        surf = (self.lam1-self.lam2)/self.lam2
        surf = surf / surf.sum()
        return np.sum(surf, axis=0)
        
    def lagprofile(self):
        surf = (self.lam1-self.lam2)/self.lam2
        surf = surf / surf.sum()
        return np.sum(surf, axis=1)
    
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
    
    # def report(self):
    #     """
    #     Report the mesurement in tabular form.
    #     """
    #     toprin
        
        
    # I/O stuff  

    def save(self,filename):
        """
        Save Measurement for future referral
        """
        io.save(self,filename)
                       
    def copy(self):
        return io.copy(self)  
            
    
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
            kwargs['vals'] = self.snrsurf
            kwargs['title'] = r'$(\lambda_1 - \lambda_2) / 2\lambda_2$'
        
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
        ax.set_ylabel(r'Fast Direction ($^\circ$)')
        ax.set_xlabel('Delay Time (' + self.units + ')')
        
        # confidence region
        if 'conf95' in kwargs and kwargs['conf95'] == True:
            ax.contour(self.lags,self.degs,self.lam2,levels=[self.lam2_95()])
            
        # marker
        if 'marker' in kwargs and kwargs['marker'] == True:
            ax.errorbar(self.lag,self.fast,xerr=self.dlag,yerr=self.dfast)

        ax.set_xlim([self.lags[0,0], self.lags[-1,0]])
        ax.set_ylim([self.degs[0,0], self.degs[0,-1]])
    
        # optional title
        if 'title' in kwargs:
            ax.set_title(kwargs['title']) 
            
        # add info in text box
        if 'info' in kwargs and kwargs['info'] == True:
            textstr = '$\phi=%.1f\pm%.1f$\n$\delta t=%.2f\pm%.2f$'%\
                        (self.fast,self.dfast,self.lag,self.dlag)
            # place a text box in upper left in axes coords
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)
                    
        return ax
        
    def plot_profiles(self,**kwargs):
        # Error analysis
        fig,ax = plt.subplots(2)
        ax0 = plt.subplot(121)
        ax1 = plt.subplot(122)

        ax0.plot(self.degs[0,:],self.fastprofile())
        ax0.axvline(self.fast)
        ax0.axvline(self.fast-2*self.dfast,alpha=0.5)
        ax0.axvline(self.fast+2*self.dfast,alpha=0.5)
        ax0.set_title('fast direction')

        ax1.plot(self.lags[:,0],self.lagprofile())
        ax1.axvline(self.lag)
        ax1.axvline(self.lag-2*self.dlag,alpha=0.5)
        ax1.axvline(self.lag+2*self.dlag,alpha=0.5)
        ax1.set_title('lag direction')

        plt.show()

    # Report
    
    class Report:        
        """
        Handle reporting of measurement.
        """
        def __init__(self):
            self.choose = []

    # Comparison
    
    def __eq__(self, other) :
        # check same class
        if self.__class__ != other.__class__: return False
        # check same keys
        if set(self.__dict__) != set(other.__dict__): return False
        # check same values
        for key in self.__dict__.keys():
            if not np.all( self.__dict__[key] == other.__dict__[key]): return False
        # if reached here then the same
        return True
        

        
