# -*- coding: utf-8 -*-
"""
The measurement class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core, core3d, io
from ..core.pair import Pair
from ..core.window import Window
# from . import eigval, rotcorr, transmin, sintens

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os.path


class Measure:
    
    """
    Base measurement class        
    """
    
    def __init__(self,*args,**kwargs):
        
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
        self.__slags = np.unique( core.time2samps( lags, self.delta, mode='even')).astype(int)

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
        self.__degs = degs
        
        # self.lags, self.degs = np.meshgrid(self.__slags * self.delta, self.__degs)
        self.degs, self.lags = np.meshgrid(self.__degs, self.__slags * self.delta)

        # receiver correction
        self.rcvcorr = None
        if ('rcvcorr' in kwargs):
            if not isinstance(kwargs['rcvcorr'],tuple): raise TypeError('rcvcorr must be tuple')
            if len(kwargs['rcvcorr']) != 2: raise Exception('rcvcorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['rcvcorr']
            samps = core.time2samps(lag, self.delta, 'even')
            self.__rcvcorr = (deg, samps)
            self.rcvcorr = kwargs['rcvcorr']

        # source correction
        self.srccorr = None
        if ('srccorr' in kwargs):
            if not isinstance(kwargs['srccorr'],tuple): raise TypeError('srccorr must be tuple')
            if len(kwargs['srccorr']) != 2: raise Exception('srccorr must be length 2')
            # convert time shift to nsamples -- must be even
            deg, lag = kwargs['srccorr']
            samps = core.time2samps(lag, self.delta, 'even')
            self.__srccorr = (deg, samps)
            self.srccorr = kwargs['srccorr']
                
    # Common methods
    
    def gridsearch(self, func, **kwargs):
        
        """
        Grid search for splitting parameters applied to data using the function defined in func
        rcvcorr = receiver correction parameters in tuple (fast,lag) 
        srccorr = source correction parameters in tuple (fast,lag) 
        """
        
        # avoid using "dots" in loops for performance
        rotate = core.rotate
        lag = core.lag
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
                
        # rotate to polaristation (needed for tranverse min)
        if 'mode' in kwargs and kwargs['mode'] == 'rotpol':
            pol = self.data.pol
            def rotpol(x, y, ang):
                # rotate to pol
                x, y = rotate(x, y, pol-ang)
                return x, y
        else:
            def rotpol(x, y, ang):
                return x, y
        
        # actual inner loop function   
        def getout(x, y, ang, shift):
            # remove shift
            x, y = lag(x, y, -shift)
            x, y = srccorr(x, y, ang)
            x, y = chop(x, y, window=self.data.window)
            x, y = rotpol(x, y, ang)
            return func(x, y)

                    
        # Do the grid search
        prerot = [ (rotate(x, y, ang), ang) for ang in self.__degs ]
        
        out = [ [ getout(data[0], data[1], ang, shift) for shift in self.__slags ]
                for (data,ang) in prerot  ]
                               
        return out
        
    def gridsearch3d(self, func, **kwargs):
        
        """
        Grid search for splitting parameters applied to data using the function defined in func
        rcvcorr = receiver correction parameters in tuple (fast,lag) 
        srccorr = source correction parameters in tuple (fast,lag) 
        """
        
        # avoid using "dots" in loops for performance
        rotate = core3d.rotate
        lag = core3d.lag
        chop = core3d.chop
        unsplit = core3d.unsplit
        
        # ensure trace1 at zero angle
        copy = self.data.copy()
        copy.rotate2ray()
        x, y, z = copy.x, copy.y, copy.z
        
        # pre-apply receiver correction
        if 'rcvcorr' in kwargs:
            rcvphi, rcvlag = self.__rcvcorr
            x, y, z = unsplit(x, y, z, rcvphi, rcvlag)
                            
        ######################                  
        # inner loop function
        ######################
    
        # source correction  
        
        if 'srccorr' in kwargs:
            srcphi, srclag = self.__srccorr
            def srccorr(x, y, z, ang):
                x, y, z = unsplit(x, y, z, srcphi-ang, srclag)
                return x, y, z
        else:
            def srccorr(x, y, z, ang):
                return x, y, z
                
        # rotate to polaristation (needed for tranverse min)
        if 'mode' in kwargs and kwargs['mode'] == 'rotpol':
            pol = self.data.pol
            def rotpol(x, y, z, ang):
                # rotate to pol
                x, y, z = rotate(x, y, z, pol-ang)
                return x, y, z
        else:
            def rotpol(x, y, z, ang):
                return x, y, z
        
        # actual inner loop function   
        def getout(x, y, z, ang, shift):
            # remove shift
            x, y, z = lag(x, y, z, -shift)
            x, y, z = srccorr(x, y, z, ang)
            x, y, z = chop(x, y, z, window=self.data.window)
            x, y, z = rotpol(x, y, z, ang)
            return func(x, y, z)
                    
        # Do the grid search
        prerot = [ (rotate(x, y, z, ang), ang) for ang in self.__degs ]
        
        out = [ [ getout(data[0], data[1], data[2], ang, shift) for shift in self.__slags ]
                for (data,ang) in prerot  ]
                               
        return out
            

            
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
        return self.data_corr().get_pol()
        
    def snr(self):
        """Restivo and Helffrich (1999) signal to noise ratio"""
        d = self.srcpoldata_corr().chop()
        return core.snrRH(d.x,d.y)
                
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
        srcpoldata.set_labels(['srcpol','trans','ray'])
        return srcpoldata
        
    def srcpoldata_corr(self):
        srcpoldata_corr = self.data_corr()        
        srcpoldata_corr.rotateto(self.srcpol())
        srcpoldata_corr.set_labels(['srcpol','trans','ray'])
        return srcpoldata_corr
        
    def fastdata(self):
        """Plot fast/slow data."""
        fastdata = self.data.copy()
        fastdata.rotateto(self.fast)
        fastdata.set_labels(['fast','slow','ray'])
        return fastdata

    def fastdata_corr(self):
        fastdata_corr = self.data_corr()
        fastdata_corr.rotateto(self.fast)
        fastdata_corr.set_labels(['fast','slow','ray'])
        return fastdata_corr
            
    # F-test utilities
    
    def ndf(self):
        """Number of degrees of freedom."""
        d = self.srcpoldata_corr().chop()
        return core.ndf(d.y)
    
    def get_errors(self,surftype=None):
        """
        Return dfast and dtlag.

        These errors correspond to one sigma in the parameter estimate.

        Calculated by taking a quarter of the width of 95% confidence region (found using F-test).
        """

        # search interval steps
        lag_step = self.lags[1,0] - self.lags[0,0]
        fast_step = self.degs[0,1] - self.degs[0,0]

        # Find nodes where we fall within the 95% confidence region
        if surftype == 'max':
            confbool = self.errsurf >= self.conf_95()
        elif surftype == 'min':
            confbool = self.errsurf <= self.conf_95()
        else:
            raise ValueError('surftype must be min or max')

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
            ax.contour(self.lags,self.degs,self.errsurf,levels=[self.conf_95()])
            
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
        

        
