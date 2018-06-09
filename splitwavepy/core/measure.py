# -*- coding: utf-8 -*-
"""
The measurement class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, io
from .data import SplitWave
# from .bootstrap import Bootstrap

# from ..core import core, core3d, io
# from ..core.pair import Pair
# from ..core.window import Window
# from . import eigval, rotcorr, transmin, sintens

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# import os.path


class Py(SplitWave):
    
    """
    Measure shearwave splitting on a SplitWave object. 
    """
    
    def __init__(self, SplitWave, **kwargs):

        # Settings
        settings = {}
        settings['plot'] = False
        settings['report'] = True
        settings['bootstrap'] = True
        settings['rcvcorr'] = None
        settings['srccorr'] = None
        
        settings.update(kwargs) # update using kwargs
        # backup keyword args
        self._settings = settings
        

        _set_degs(**settings)
        _set_lags(**settings)
        self._rcvcorr = settings['rcvcorr']
        self._srccorr = settings['srccorr']
        


            
        # Grid Search
        self.covmap = self._gridcov(**kwargs)
        
        # Inspect the covariance map
        # Silver and Chan
        self.lam1, self.lam2 = self._silver_and_chan(**kwargs)
        sc_loc = core.max_idx(self.lam1 / self.lam2)
        self.sc_fast, self.sc_lag = self.degmap[sc_loc], self.lagmap[sc_loc]

        # Cross-correlation
        self.xc = core.covmap2rho()
        xc_loc = core.max_idx(self.xc)
        self.xc_fast, self.xc_lag =  self.degmap[xc_loc], self.lagmap[xc_loc]
        
        # Splitting Intensity
        self.splintensity = self._splitting_intensity()
        
        # Name
        self.name = 'Untitled'
        if 'name' in kwargs: self.name = kwargs['name']
            

        
        # Implement Settings
        if 'bootstrap' == True : m.bootstrap(**kwargs)            
        if 'report'    == True : m.report(**kwargs)
        if 'plot'      == True : m.plot(**kwargs)

    #===================
    # Special Properties
    #===================
    
    #
    # # _data
    # @property
    # def _data(self):
    #     return self.__data
    #

    # _degs
    @property
    def _degs(self):
        return self.__degs
        
    @_degs.setter
    def _degs(self, degs):
        if not isinstance(degs, np.ndarray):
            raise TypeError('degs must be numpy array')
        if degs.ndim != 1:
            raise ValueError('degs must be 1-d')
        self.__degs = degs
        self.__rads = np.radians(degs)
        
    def _set_degs(self, **kwargs):
        """return numpy array of degs to explore"""
        settings = {}
        settings['mindeg'] = -90
        settings['maxdeg'] = 90
        settings['ndegs']  = 90
        settings.update(kwargs)
        self._degs = np.linspace(settings['mindeg'],
                                 settings['maxdeg'],
                                 settings['ndegs'],
                                 endpoint=False)
        
    #_lags
    @property
    def _lags(self):
        return self.__lags
        
    @_lags.setter
    def _lags(self, lags):
        if not isinstance(lags, np.ndarray):
            raise TypeError('lags must be numpy array')
        if lags.ndim != 1:
            raise ValueError('lags must be 1-d')
        self.__slags = np.unique(core.time2samps(lags, self._delta, mode='even')).astype(int)
        self.__lags = self.__slags * self._delta
    
    def _set_lags(self, **kwargs):
        """return numpy array of lags to explore"""
        settings = {}
        settings['minlag'] = 0
        settings['maxlag'] = self.wwidth() / 4
        settings['nlags']  = 40
        settings.update(kwargs)
        self._lags = np.linspace(settings['minlag'],
                                 settings['maxlag'],
                                 settings['nlags'],
                                 endpoint = True)
                                 
    #_grid
    @property
    def _grid(self):
        return np.meshgrid(self._degs, self._lags)
                   
    #_rcvcorr    
    @property
    def _rcvcorr(self):
        return self.__rcvcorr
        
    @_rcvcorr.setter
    def _rcvcorr(self, rcvcorr):
        if rcvcorr == None:
            self.__rcvslag = None
            self.__rcvcorr = None
        if isinstance(rcvcorr, tuple) and len(rcvcorr) == 2: 
            deg, lag = rcvcorr
            slag = core.time2samps(lag, self._delta, 'even')
            self.__rcvslag = (deg, slag)
            self.__rcvcorr = (deg, slag * self._delta)
        else: raise ValueError('rcvcorr not understood.')
                
    #_srccorr
    @property
    def _srccorr(self):
        return self.__srccorr
        
    @_srccorr.setter
    def _srccorr(self, srccorr):
        if srccorr == None:
            self.__rcvslag = None
            self.__srccorr = None
        if isinstance(srccorr, tuple) and len(srccorr) == 2: 
            deg, lag = srccorr
            slag = core.time2samps(lag, self._delta, 'even')
            self.__rcvslag = (deg, slag)
            self.__srccorr = (deg, slag * self._delta)
        else: raise ValueError('srccorr not understood.')      

                

    
                
    # Common methods
            
    def _gridcov(self, **kwargs):       
        """
        Grid search for splitting parameters applied to self.SplitWave using the function defined in func
        rcvcorr = receiver correction parameters in tuple (fast,lag) 
        srccorr = source correction parameters in tuple (fast,lag) 
        """
        
        x, y = np.copy(self.x), np.copy(self.y)
        w0, w1 = self._w0(), self._w1()
                
        # receiver correction
        if 'rcvcorr' in kwargs:
            rcvphi, rcvlag = self.__rcvcorr
            x, y = core.unsplit(x, y, rcvphi, rcvlag)

        # source correction
        if 'srccorr' in kwargs:
            raise Exception('Not implemented.')
            # srcphi, srclag = self.__srccorr
            # return core.gridcov_srcorr(x, y, w0, w1, degs, slags, srcphi, srclag)
        
        return core.gridcov(x, y, w0, w1, self._degs, self._slags)
        
    def _silver_and_chan(self):
        if 'pol' in kwargs:
            raise Exception('Not implemented.')
            #lam1, lam2 = core.covmap2polvar(self.covmap, pol)
        else:
            # use eigen analysis
            lam1, lam2 = core.covmap2eigvals(self.covmap)
        return lam1, lam2
        
    def _correlation(self):
        return core.covmap2rho(self.covmap)
        
    # utility
    



                    
    # METHODS 
    #--------    
    
    def report(self):
        raise Exception('Not implemented.')
    
    def srcpol(self):
        # recover source polarisation
        if 'pol' in self.kwargs:
            return self.kwargs['pol']
        else:
            return self.SplitWave_corr().estimate_pol()
        
    def snr(self):
        """Restivo and Helffrich (1999) signal to noise ratio"""
        d = self.srcpoldata_corr()._chop()
        return core.snrRH(d.x, d.y)
                
    # data views
    
    def data_corr(self):        
        # copy data     
        data_corr = self.copy()
        # rcv side correction     
        if self.rcvcorr is not None:
            data_corr = data_corr.unsplit(*self.rcvcorr)    
        # target layer correction
        data_corr = data_corr.unsplit(self.fast, self.lag)  
        # src side correction
        if self.srccorr is not None:
            data_corr = data_corr.unsplit(*self.srccorr)
        return data_corr

    def srcpoldata(self):
        srcpoldata = self.copy()
        srcpoldata._rotateto(self.srcpol())
        srcpoldata.set_labels(['srcpol', 'trans', 'ray'])
        return srcpoldata
        
    def srcpoldata_corr(self):
        srcpoldata_corr = self.SplitWave_corr()        
        srcpoldata_corr._rotateto(self.srcpol())
        srcpoldata_corr.set_labels(['srcpol', 'trans', 'ray'])
        return srcpoldata_corr
        
    def fastdata(self):
        """Plot fast/slow data."""
        fastdata = self.copy()
        fastdata._rotateto(self.fast)
        fastdata.set_labels(['fast', 'slow', 'ray'])
        return fastdata

    def fastdata_corr(self):
        fastdata_corr = self.SplitWave_corr()
        fastdata_corr._rotateto(self.fast)
        fastdata_corr.set_labels(['fast', 'slow', 'ray'])
        return fastdata_corr
            
    # F-test utilities
    
    def ndf(self):
        """Number of degrees of freedom."""
        x, y = self.srcpoldata_corr()._chopdata()
        return core.ndf(y)
    
    def get_errors(self, surftype=None):
        """
        Return dfast and dtlag.

        These errors correspond to one sigma in the parameter estimate.

        Calculated by taking a quarter of the width of 95% confidence region (found using F-test).
        """

        # search interval steps
        lag_step = self.lags[1] - self.lags[0]
        fast_step = self.degs[1] - self.degs[0]

        # Find nodes where we fall within the 95% confidence region
        
        if surftype == 'max':
            confbool = self.errsurf >= self.conf95level
        elif surftype == 'min':
            confbool = self.errsurf <= self.conf95level
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
        cyclic = np.hstack((fastbool, fastbool))
        lengthFalse = np.diff(np.where(cyclic)).max() - 1
        # shortest line that contains ALL true values is then:
        lengthTrue = fastbool.size - lengthFalse
        fdfast = lengthTrue * fast_step * 0.25

        # return
        return fdfast, fdlag 
        
    # splitting intensity
    def splitting_intensity(self, pol=None, **kwargs):
        """
        Calculate the splitting intensity as defined by Chevrot (2000).
        """        
        
        # settings
        
        if 'pol' not in kwargs:
            raise Exception('pol must be specified')            
        copy = self.copy()
        copy._rotateto(pol)
        copy.x = np.gradient(copy.x)
        rdiff, trans = copy._chopdata()
        s = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
        return s
        
    # bootstrap utilities
    # to do: implement block bootstrapping for time series data.
    
    def _bootstrap_samp(self, x, y):
        """Calculate a single bootstrap statistic on one resampling of the x, y data."""
        return self._bootstrap_stat(*core.bootstrap_resamp(x, y))
    
    def _bootstrap_loop(self, *args, n=5000, **kwargs):
        """Calculate many bootstrap statistics on n resamplings of the data."""
        # ensure data prepped (in correct orientation and windowed) appropriately
        if len(args) == 2:
            x, y = args[0], args[1]
        else:
            x, y = self._bootstrap_prep()
        # calculate bootstrap values
        bootstrap_vals = np.asarray([ self._bootstrap_samp(x, y) for ii in range(n) ])
        return bootstrap_vals
        
    def _bootstrap_grid(self, **kwargs):
        return self.gridsearch(self._bootstrap_loop, **kwargs)
        
    def estimate_pdf(self, polar=False, **kwargs):
        x, y = self._bootstrap_prep()
        vals = self._bootstrap_loop(x, y, **kwargs)
        kde = core.kde(vals)
        ravmap = np.ravel(self.vals())
        pdf = kde.pdf(ravmap).reshape(self.vals().shape)
        
        # if polar == True:
        #     # geometric correction to account for variation
        #     # in sampling density if grid is converted to polar coordinates
        #     # (more research needed to establish whether this is correct)
        #     ###
        #     # circle area is pi*r^2
        #     # annulus area is pi*r1^2 - pi*r0^2
        #     # cell area is annulus area / number of cells
        #     # cell area is K*(r1^2-r0^2), where K is pi/number of cells.
        #     # K is just a constant scaling factor that can be removed at the end.
        #     aa = np.arange(self.lags.size)
        #     bb = (aa+1)**2 - aa**2
        #     _, cc = np.meshgrid(self.degs, bb)
        #     pdf = pdf * cc
            
        # normalise so that whole surface weighs 1
        pdf = pdf / np.sum(pdf)
        
        return pdf
        
    def _pdf_conf95(self, pdf):
        return core.val_at_alpha(pdf.flatten(),0.05)
        
    # def pdf(self, **kwargs):
    #     return self.estimate_pdf(**kwargs)

    # error propagating corrections

    def _correction_variance(self, rcvinfo=None, srcinfo=None, n=100, m=100):
        """Propagate errors in receiver and/or source correction.
        rcvinfo = (fast, dfast, lag, dlag)
        srcinfo = (fast, dfast, lag, dlag)
        n = number of trial corrections
        m = number of bootstrap subsamples per trial

        Uses bootstrapping to calculate 95% confidence level.
        Trials receiver corrections randomly drawn from a normal distribution."""

        def _get_data(rcvcorr=None, srccorr=None):
            """same as data_corr but user can change corrections"""
            # copy data
            data_corr = self.copy()
            # rcv side correction
            if rcvcorr is not None:
                data_corr = data_corr.unsplit(*rcvcorr)
            # target layer correction
            data_corr = data_corr.unsplit(self.fast, self.lag)
            # src side correction
            if srccorr is not None:
                data_corr = data_corr.unsplit(*srccorr)                
            # ensure orientation of data is appropriate for func
            if self.func == core.transenergy:
                return data_corr._chopdata().rotatetto(self.srcpol())
            elif (self.func == core.crosscorr) or (self.func == core.pearson):
                return data_corr._chopdata()._rotateto(self.fast)
            else:
                return data_corr._chopdata()

        def _draw_corr(info):
            """draw a correction from a gaussian distribution"""
            if info is None: return None
            fast = np.random.normal(info[0], info[1])
            lag = np.random.normal(info[2], info[3])
            return fast, lag

        datafeed = ( _get_data(rcvcorr=_draw_corr(rcvinfo), 
                               srccorr=_draw_corr(srcinfo))
                               for ii in range(n) )
                               
        vals =  [ [ self._bootstrap_samp(*data) for ii in range(m) ] 
                    for data in datafeed ]
                    
        return np.asarray(vals).flatten()
        
        #
        # # vals = np.asarray([ self.func(*data) for data in datafeed ])
        #
        # if self.func == core.transenergy:
        #     # get minimimum energy and keep upper half
        #     return np.sort(vals[:,1])[int(m/2):-1]
        # elif (self.func == core.crosscorr) or (self.func == core.pearson):
        #     # get coefficients and keep lower half
        #     return np.sort(vals)[0:int(m/2)]
        # elif self.func == core.eigvalcov:
        #     # get minimum eigenvalue and keep upper half
        #     return np.sort(vals[:,0])[int(m/2):-1]

        
        # def _bootstrap(data):
        #     """Bootstrap the data after trial correction applied"""
        #     # keep only lower/upper half to make one-sided distribution
        #     vals = np.asarray([ self.func(*core.bootstrap_resamp(*data)) for ii in range(m) ])
        #     if self.func == core.transenergy:
        #         # get minimimum energy and keep upper half
        #         return np.sort(vals[:,1])[int(m/2):-1]
        #     elif (self.func == core.crosscorr) or (self.func == core.pearson):
        #         # get coefficients and keep lower half
        #         return np.sort(vals)[0:int(m/2)]
        #     elif self.func == core.eigvalcov:
        #         # get minimum eigenvalue and keep upper half
        #         return np.sort(vals[:,0])[int(m/2):-1]
        #
        # return [ _bootstrap(data) for data in datafeed ]
            
        

            
            
    # bootstrap a la Sandvol and Hearn
            
    def _renoise(self, **kwargs):
        """
        Return data with new noise sequence
        """
        # copy original, corrected, data
        newself = self.copy()
        bs = self.SplitWave_corr()
        origang = bs.cmpangs()[0]
        # replace noise sequence
        bs._rotateto(self.srcpol())
        bs.y = core.resample_noise(bs.y)
        bs._rotateto(origang)
        # reapply splitting
        # src side correction
        if self.srccorr is not None: bs = bs._split(*self.srccorr)
        # target layer correction
        bs = bs._split(self.fast, self.lag)
        # rcv side correction
        if self.rcvcorr is not None: bs = bs._split(*self.rcvcorr)
        newself.SplitWave = bs
        return newself

    def _bootstrap_sandhgrid(self, **kwargs):

        return ( self._renoise(**kwargs) for x in range(kwargs['n']) )
        # return ( newself.gridsearch(**kwargs) for newself in newselffeed )
        #

    # def bootstrap(self, **kwargs):
    #     return Bootstrap(self)
        
    # "squashed" profiles
    
    def fastprofile(self, **kwargs):
        if 'vals' not in kwargs:
            raise Exception('vals must be specified')
        surf = kwargs['vals']
        surf = surf / surf.sum()
        return np.sum(surf, axis=0)
        
    def lagprofile(self, **kwargs):
        if 'vals' not in kwargs:
            raise Exception('vals must be specified')
        surf = kwargs['vals']
        surf = surf / surf.sum()
        return np.sum(surf, axis=1)
        


    
    # Output
    
    # def report(self):
    #     """
    #     Report the mesurement in tabular form.
    #     """
    #     toprin
        
        
    # I/O stuff

    def save(self, filename):
        """
        Save Pyment for future referral
        """
        io.save(self, filename)
        
    def load(self, filename):
        return io.load(self, filename)

    def copy(self):
        return io.copy(self)
        
    # spit out the answer
        
    def report(self, **kwargs):
        """Prints fast, lag, dfast, dlag to screen/stdout."""
        print('fast'.rjust(10), 'dfast'.rjust(9), 'lag'.rjust(9), 'dlag'.rjust(9))
        print('{0:10.2f}{1:10.2f}{2:10.2f}{3:10.2f}'.format(self.fast, self.dfast, self.lag, self.dlag))
            
    
    # Plotting
    
    def _plot(self, **kwargs):
        
        if 'vals' not in kwargs:
            raise Exception('vals must be specified')
            
        # if kwargs['vals'] == 'pdf':
        #     kwargs['vals'] = self.estimate_pdf()
        #     kwargs['title'] = r'Probability Density'
          
        # setup figure and subplots
        fig = plt.figure(figsize=(12,6)) 
        
        gs = gridspec.GridSpec(3, 3,
                           width_ratios=[2,1,3]
                           )
        ax0 = plt.subplot(gs[0,0:2])                     
        ax1 = plt.subplot(gs[1,0])
        ax2 = plt.subplot(gs[1,1])
        ax3 = plt.subplot(gs[2,0])
        ax4 = plt.subplot(gs[2,1])
        ax5 = plt.subplot(gs[:,2])

        orig = self.srcpoldata()._chop()
        corr = self.srcpoldata_corr()._chop()
                
        # get axis scaling
        lim = np.abs(corr.data()).max() * 1.1
        ylim = [-lim, lim]
        
        # long window data
        self._ptr(ax0, ylim=ylim, **kwargs)

        # original
        orig._ptr(ax1, ylim=ylim, **kwargs)
        orig._ppm(ax2, lims=ylim, **kwargs)
        
        # corrected
        corr._ptr(ax3, ylim=ylim, **kwargs)
        corr._ppm(ax4, lims=ylim, **kwargs)
        
        # add marker and info box by default
        if 'marker' not in kwargs: kwargs['marker'] = True
        if 'info' not in kwargs: kwargs['info'] = True
        if 'conf95' not in kwargs: kwargs['conf95'] = True
        self._psurf(ax5,**kwargs)
        
        # title
        if 'name' in kwargs:
            plt.suptitle(kwargs['name'])
                    
        # neaten
        plt.tight_layout()
        
        # save or show
        if 'file' in kwargs:
            plt.savefig(kwargs['file'])
        else:
            plt.show()
    
        
    def _psurf(self, ax, **kwargs):
        """
        Plot an error surface.
    
        **kwargs
        - cmap = 'magma'
        - vals = (M.lam1-M.lam2) / M.lam2
        """
    
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'magma'
    
        if 'vals' not in kwargs:
            raise Exception('vals must be specified')
            
        # error surface
        cax = ax.contourf(self.lagmap, self.degmap, kwargs['vals'], 26, cmap=kwargs['cmap'])
        cbar = plt.colorbar(cax)
        ax.set_ylabel(r'Fast Direction ($^\circ$)')
        ax.set_xlabel('Delay Time (' + self.units + ')')
        
        # confidence region
        # if 'conf95' in kwargs and kwargs['conf95'] == True:
        #     ax.contour(self.lagmap, self.degmap, self.errsurf, levels=[self.conf95level],
        #             colors='r', alpha=.5, linewidths=3)
            
        # marker
        if 'marker' in kwargs and kwargs['marker'] == True:
            ax.errorbar(self.lag, self.fast, xerr=self.dlag, yerr=self.dfast)

        ax.set_xlim([laggrid[0,0], laggrid[-1,0]])
        ax.set_ylim([deggrid[0,0], deggrid[0,-1]])
    
        # optional title
        if 'title' in kwargs:
            ax.set_title(kwargs['title']) 
            
        # add info in text box
        if 'info' in kwargs and kwargs['info'] == True:
            textstr = '$\phi=%.1f\pm%.1f$\n$\delta t=%.2f\pm%.2f$'%\
                        (self.fast, self.dfast, self.lag, self.dlag)
            # place a text box in upper left in axes coords
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                    verticalalignment='top', bbox=props)
        
        if 'ppm' in kwargs and kwargs['ppm'] is True:
            sublags = self.lags[0:-1:int(self.lags.size/6)]
            subdegs = self.degs[0:-1:int(self.degs.size/6)]
            sublags = sublags + (self.lags[-1]-sublags[-1]) / 2
            subdegs = subdegs + (self.degs[-1]-subdegs[-1]) / 2
            x, y = self.SplitWave_corr()._chopdata()   
            lagtot = self.lags[-1] - self.lags[0]
            degtot = self.degs[-1] - self.degs[0]
            boost = 10 * lagtot / np.max((x**2 + y**2)**.5)      
            for fast in subdegs:
                for lag in sublags:
                    x, y = self.unsplit(fast, lag)._chopdata()
                    ax.plot(lag + y*boost/degtot, fast + x*boost/lagtot, color='w',alpha=0.5)

                    
        return ax
        

        
    def plot_profiles(self,**kwargs):
        # Error analysis
        fig,ax = plt.subplots(2)
        ax0 = plt.subplot(121)
        ax1 = plt.subplot(122)

        ax0.plot(self.degs[0,:], self.fastprofile())
        ax0.axvline(self.fast)
        ax0.axvline(self.fast-2*self.dfast, alpha=0.5)
        ax0.axvline(self.fast+2*self.dfast, alpha=0.5)
        ax0.set_title('fast direction')

        ax1.plot(self.lags[:,0], self.lagprofile())
        ax1.axvline(self.lag)
        ax1.axvline(self.lag-2*self.dlag, alpha=0.5)
        ax1.axvline(self.lag+2*self.dlag, alpha=0.5)
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
        

        
