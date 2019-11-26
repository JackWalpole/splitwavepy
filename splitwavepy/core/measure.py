# -*- coding: utf-8 -*-
"""
The measurement class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, io
from .data import Data
# from .bootstrap import Bootstrap

# from ..core import core, core3d, io
# from ..core.pair import Pair
# from ..core.window import Window
# from . import eigval, rotcorr, transmin, sintens

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats


# import os.path


class Meas(Data):
    
    """
    Measure shearwave splitting on a Data object.
    
    Usage: Meas(Data, **options)
    
    Calculates the covariance matrix of unsplit data to search for splitting
    parameters that best remove the splitting.
    
    Splitting is characterised as either:
    
    -  the ratio of the energy in the polarisation direction
    to the energy in the transverse direction (Silver and Chan, 1991).
    
    - the correlation of the fast and slow wavelets (e.g. Ando, XXXX).
    
    
    
    Uses the Fourier transform to do the grid search by default, unless the 
    
    **options
    ---------
    
    maxlag  
    ndegs
    
    Methods:
    
    plot
    
    
    
    """
    
    def __init__(self, Data, **kwargs):
        
        # The Data object
        self.__data = Data.rotateto(0)

        # Default Settings (These can be overidden using keywords)
        settings = {}
        settings['plot'] = False
        settings['report'] = False
        settings['bootstrap'] = True 
        settings['rcvcorr'] = None
        settings['srccorr'] = None
        settings['taper'] = 0.2
        settings['fft'] = True
        
        # degs settings
        settings['pol'] = None
        settings['ndegs'] = 180
        settings['maxlag'] = None
        settings.update(kwargs) # override defaults using kwargs
        self._settings = settings # backup settings
        
        # Implement setting
        # self._set_degs(**settings)
        # self._set_lags(**settings)
        self._pol = settings['pol']
        self._rcvcorr = settings['rcvcorr']
        self._srccorr = settings['srccorr']
        self._ndegs = settings['ndegs']
        self._maxlag =  settings['maxlag']
        self._taper = settings['taper']
        self._bootstrap = settings['bootstrap']
            
        # Grid Search
        # self._covmap = self._gridcov()
        if settings['fft']: self._covmap = self._gridcovfreq()
        else : self._covmap = self.gridcov() # logic here to reduce number of lag steps
        
        
        # self.sc = self.silver_chan()
        # self.xc = self.correlation()
        # self.q = core.q(self.sc.fast, self.sc.lag, self.xc.fast, self.xc.lag)

        # # Splitting Intensity
        # self.splintensity = self._splitting_intensity()
        #
        # # Name
        # self.name = 'Untitled'
        # if 'name' in kwargs: self.name = kwargs['name']
            
        # print(self.sc.srcpol)

        # Implement Settings
        # if 'bootstrap' == True : m.bootstrap(**kwargs)
        if settings['report'] == True : self.report(**kwargs)
        # if 'plot'      == True : m.plot(**kwargs)

    #===================
    # Special Properties
    #===================
    

    # _data
    @property
    def data(self):
        return self.__data
        
    def data_corr(self, fast, lag):
        return self.data.unsplit(fast, lag)
        
    def srcpoldata(self, pol=None):
        if pol == None: pol = self.data._estimate_pol()
        srcpoldata = self.data.rotateto(pol)
        srcpoldata._set_labels(['pol', 'trans'])
        return srcpoldata

    def srcpoldata_corr(self, fast, lag, pol=None):
        if pol == None: pol = self.data._estimate_pol()
        srcpoldata_corr = self.data_corr(fast, lag).rotateto(pol)
        srcpoldata_corr._set_labels(['pol', 'trans'])
        return srcpoldata_corr

    def fastdata(self, fast):
        """Plot fast/slow data."""
        fastdata = self.data.rotateto(fast)
        fastdata._set_labels(['fast', 'slow'])
        return fastdata

    def fastdata_corr(self, fast, lag):
        fastdata_corr = self.data_corr(fast, lag).rotateto(fast)
        fastdata_corr._set_labels(['fast', 'slow'])
        return fastdata_corr
        
    #_ndegs
    @property
    def _ndegs(self):
        return self.__ndegs
        
    @_ndegs.setter
    def _ndegs(self, ndegs):
        self.__ndegs  = int(ndegs)
        
    @property
    def _degs(self):
        return np.linspace(0, 180, self._ndegs, endpoint=False)
        
    #_maxlag
    @property
    def _maxlag(self):
        return self.__maxlag
        
    @_maxlag.setter
    def _maxlag(self, maxlag):
        if maxlag is None:
            maxlag = self.data.wwidth() / 4
        maxslag = core.time2samps(maxlag, self.data._delta)
        maxlag = maxslag * self.data._delta
        self.__nslags = maxslag + 1
        self.__maxlag = maxlag
        
    @property
    def _nlags(self):
        return self.__nslags

    @property
    def _slags(self):
        return np.arange(self.__nslags)
    
    @property
    def _lags(self):
        return self._slags * self.data._delta



    @property
    def _pol(self):
        return self.__pol

    @_pol.setter
    def _pol(self, pol):
        self.__pol = pol
        
    # @_pol.setter
    # def _pol(self, pol=None):
    #     if pol is None: self.__pol = self._guesspol()
    #     else: self.__pol = float(pol)
        

    #_degs
    # @property
    # def _degs(self):
    #     return np.linspace(0, 180, settings['ndegs'],
    #                              endpoint=False)

    # @_degs.setter
    # def _degs(self, degs):
    #     if isinstance(degs, np.ndarray) and degs.ndim == 1:
    #         self.__degs = degs
    #         # self.__rads = np.radians(degs)
    #     else: raise ValueError('degs not understood.')

    # def _set_degs(self, **kwargs):
    #     """return numpy array of degs to explore"""
    #     settings = {}
    #     settings['mindeg'] = -90
    #     settings['maxdeg'] = 90
    #     settings['ndegs']  = 90
    #     settings.update(kwargs)
    #     self._degs = np.linspace(settings['mindeg'],
    #                              settings['maxdeg'],
    #                              settings['ndegs'],
    #                              endpoint=False)
        
    #_lags
    # @property
    # def _lags(self):
    #     return self.__lags
    #
    # @_lags.setter
    # def _lags(self, lags):
    #     if isinstance(lags, np.ndarray) and lags.ndim == 1:
    #         self.__slags = np.unique(core.time2samps(lags, self.data._delta, mode='even')).astype(int)
    #         self.__lags = self.__slags * self.data._delta
    #     else: raise ValueError('lags not understood.')
    #
    # def _set_lags(self, **kwargs):
    #     """return numpy array of lags to explore"""
    #     settings = {}
    #     settings['minlag'] = 0
    #     settings['maxlag'] = self.data.wwidth() / 4
    #     settings['nlags']  = 40
    #     settings.update(kwargs)
    #     self._lags = np.linspace(settings['minlag'],
    #                              settings['maxlag'],
    #                              settings['nlags'],
    #                              endpoint = True)
                                 

                                 
    #_grid
    @property
    def _grid(self):
        return np.meshgrid(self._lags, self._degs, indexing='ij')
                   
    #_rcvcorr    
    @property
    def _rcvcorr(self):
        return self.__rcvcorr
        
    @_rcvcorr.setter
    def _rcvcorr(self, rcvcorr):
        if rcvcorr == None:
            self.__rcvslag = None
            self.__rcvcorr = None
        elif isinstance(rcvcorr, tuple) and len(rcvcorr) == 2: 
            deg, lag = rcvcorr
            slag = core.time2samps(lag, self.data._delta, 'even')
            self.__rcvslag = (deg, slag)
            self.__rcvcorr = (deg, slag * self.data._delta)
        else: raise TypeError('rcvcorr not understood.')
                
    #_srccorr
    @property
    def _srccorr(self):
        return self.__srccorr
        
    @_srccorr.setter
    def _srccorr(self, srccorr):
        if srccorr == None:
            self.__rcvslag = None
            self.__srccorr = None
        elif isinstance(srccorr, tuple) and len(srccorr) == 2: 
            deg, lag = srccorr
            slag = core.time2samps(lag, self.data._delta, 'even')
            self.__rcvslag = (deg, slag)
            self.__srccorr = (deg, slag * self.data._delta)
        else: raise TypeError('srccorr not understood.')
        
     
        
    def splitting_intensity(self, pol=None, **kwargs):
        """
        Calculate the splitting intensity as defined by Chevrot (2000).
        """                
        # need to rotate data to source polarisation direction.
        if pol == None: pol = self.data._estimate_pol()       
        copy = self.data.rotateto(pol)
        copy.__x = np.gradient(copy.x)
        rdiff, trans = copy._chopxy()
        s = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
        return s      
    
    def silver_chan(self, alpha=0.05, pol=None):
        """Return splitting parameters plus error bars at alpha."""
        if pol == None:
            return self.ferror_min(self.lam2, alpha)
        else:
            self.pol = pol
            return self.ferror_min(self.rad2, alpha)
     
    def cross_corr(self, alpha=0.05):
         return self.ferror_max(self.zrho, alpha)
         
    def q(self):
        scfast, _, sclag, _ = self.silver_chan()
        xcfast, _, xclag, _ = self.cross_corr()
        return core.q(scfast, sclag, xcfast, xclag)
         
    # Visible methods
    
    def report(self, **kwargs):
        # scfast, scdfast, sclag, scdlag = self.silver_chan()
        # xcfast, xcdfast, xclag, xcdlag = self.cross_corr()
        print('SCfast'.rjust(10), 'SCdfast'.rjust(9), 'SClag'.rjust(9), 'SCdlag'.rjust(9),
              'XCfast'.rjust(9), 'XCdfast'.rjust(9), 'XClag'.rjust(9), 'XCdlag'.rjust(9), 
              'Q'.rjust(9), 'SI'.rjust(9))
        print('{0:10.2f}{1:10.2f}{2:10.2f}{3:10.2f}{4:10.2f}{5:10.2f}{6:10.2f}{7:10.2f}{8:10.2f}{9:10.5f}'.format(
               *self.silver_chan(), 
               *self.cross_corr(), 
               self.q(), self.splitting_intensity(**kwargs)))


    # interrogate the covmap
        
    # eigenvalues
    @property
    def _eigvals(self):
        """Silver and Chan's eigenvalues."""
        return core.covmap2eigvals(self._covmap)

    @property
    def lam1(self):
        """Big eigenvalue."""
        return self._eigvals[0]
    
    @property
    def lam2(self):
        """Small eigenvalue."""
        return self._eigvals[1]
    
    @property
    def lamrat(self):
        """Ratio of Big/Small eigenvalues."""
        lam1, lam2 = self._eigvals
        return lam1/lam2
        
    # fixed polarisation
    @property
    def _polenergy(self):
      """Radial and Transverse energy, using self._pol."""
      covrot = core.cov_rotate(self._covmap, self._pol)
      rad, trans = covrot[:,:,0,0], covrot[:,:,1,1]
      return rad, trans
      
    @property
    def rad1(self):
         """Radial energy."""
         return self._polenergy[0]
     
    @property
    def rad2(self):
         """Transverse energy."""
         return self.polenergy[1]
         
    @property
    def radrat(self):
         rad1, rad2 = self._polenergy
         return rad1/rad2
         
    @property
    def rho(self): 
        """Pearson correlation coefficient."""
        return np.abs(core.covmap2rho(self._covmap))
    
    @property
    def zrho(self):
        """Fisher transform rho."""
        return np.arctanh(self.rho)
    

        
    

    
                
    # Common methods
    
    def _fast_lag_maxloc(self, vals):
        idx = core.max_idx(vals)
        ll, dd = self._grid
        return dd[idx], ll[idx]

    def _fast_lag_minloc(self, vals):
        idx = core.min_idx(vals)
        ll, dd = self._grid
        return dd[idx], ll[idx]
    
    # def _guesspol(self):
    #     return self.sc.srcpol
            
    def _gridcov(self):       
        """
        Grid search for splitting parameters applied to self.SplitWave using the function defined in func
        rcvcorr = receiver correction parameters in tuple (fast,lag) 
        srccorr = source correction parameters in tuple (fast,lag) 
        """
        
        # ensure data rotated to zero
        data = self.data.rotateto(0)
        x, y = data.x, data.y
        w0, w1 = data._w0(), data._w1()
                
        # receiver correction
        if self._rcvcorr is not None:
            fast, slag = self._rcvslag
            x, y = core.unsplit(x, y, fast, slag)

        # source correction
        if self._srccorr is not None:
            raise NotImplementedError('Not implemented.')
            # srcphi, srclag = self.__srccorr
            # return core.gridcov_srcorr(x, y, w0, w1, degs, slags, srcphi, srclag)
        
        cov = core.gridcov(x, y, w0, w1, self.__degs, self.__slags)
        cov = core.covfreq_reshape(cov)  
        return cov
        
    def _gridcovfreq(self, **kwargs):       
        """
        Fast shear wave splitting using Fourier transform.
        rcvcorr = receiver correction parameters in tuple (fast,lag) 
        """
        
        # ensure data rotated to zero
        data = self.data.rotateto(0)
        x, y = data.x, data.y
                
        # receiver correction
        if self._rcvcorr is not None:
            fast, slag = self._rcvslag
            x, y = core.unsplit(x, y, fast, slag)

        # source correction
        if self._srccorr is not None:
            raise NotImplementedError('Not implemented.')
            # srcphi, srclag = self.__srccorr
            # return core.gridcov_srcorr(x, y, w0, w1, degs, slags, srcphi, srclag)
            
        # chop to window
        w0, w1 = data._w0(), data._w1()
        x, y = core.chop(x, w0, w1), core.chop(y, w0, w1)
        
        # taper because samples wrap around inside window and therefore
        # prudent to reduce the potential influence of samples near edge of window.
        # if 'taper' not in kwargs: kwargs['taper'] = 0.8
        x = core.taper(x, alpha=self._taper)
        y = core.taper(y, alpha=self._taper) 
        
        cov = core.gridcovfreq(x, y, ndegs=self.__ndegs, nslags=self.__nslags)
        cov = core.cov_reshape(cov)
        return cov
        
    # PDF
    
    # def likeli(vals, stat):
    #     return stat.pdf(vals)
    #
    # def pdf(vals, stat):
    #
    #     return vals / np.sum(vals)
    
    @property    
    def bslamratpdf(self, **kwargs):
        vals = self.lamrat
        fast, lag = self._fast_lag_maxloc(vals)
        bsvals = self._bslamratnodefft(fast, lag, **kwargs)
        return self._bsnormpdf(bsvals, vals)

    @property
    def bslam2pdf(self, **kwargs):
        vals = self.lam2
        fast, lag = self._fast_lag_minloc(vals)
        bsvals = self._bslam2nodefft(fast, lag, **kwargs)
        return self._bsnormpdf(bsvals, vals)
    
    @property
    def bszrhopdf(self, **kwargs):
        vals = np.arctanh(self.rho)
        fast, lag = self._fast_lag_maxloc(vals)
        bsvals = np.arctanh(self._bsrhonodefft(fast, lag, **kwargs))
        return self._bsnormpdf(bsvals, vals)
        
    def _bsnormpdf(self, bsvals, vals):
        norm = core.norm(bsvals)
        likelihood = norm.pdf(vals)
        pdf = likelihood / np.sum(likelihood)
        return pdf

    def _bskdepdf(self, bsvals, vals):
        kde = core.kde(bsvals)
        likelihood = kde.pdf(vals.flatten()).reshape((vals.shape))
        pdf = likelihood / np.sum(likelihood)
        return pdf
        
    # Bootstrap surfaces
        
    def _bscovnodefft(self, fast, lag, **kwargs):
        """bootstrap covariance matrices at a node."""
        
        if self._rcvcorr is not None:
            raise NotImplementedError('Not yet implemented.')
            
        if self._srccorr is not None:
            raise NotImplementedError('Not yet implemented.')
            
        # prepare data by backing off the splitting parameters
        # uses wraparound in time window to replicate behaviour in frequency domain.
        bsindata = self.data._wrap_unsplit(fast, lag, taper=self._taper, **kwargs)
        x, y = bsindata.x, bsindata.y
        return core.bootcov(x, y, **kwargs)
        
    # def _bscovnode(self, fast, lag, **kwargs):
    #     return self._bootcov(fast, lag)
        
    def _bslam2nodefft(self, fast, lag, **kwargs):
        return core.bscov2lam2(self._bscovnodefft(fast, lag, **kwargs))
        
    def _bslamratnodefft(self, fast, lag, **kwargs):
        return core.bscov2lamrat(self._bscovnodefft(fast, lag, **kwargs))
        
    def _bsrhonodefft(self, fast, lag, **kwargs):
        return core.bscov2rho(self._bscovnodefft(fast, lag, **kwargs))
    

    

     

        
    # F-test PDF surfaces
    
    def ndf(self, fast, lag, pol=None):
        # fast, lag = self._fast_lag_maxloc(vals)
        unsplit = self.data._wrap_unsplit_rotate_back(fast, lag)
        if pol == None: pol = unsplit._estimate_pol()
        x, y = unsplit.rotateto(pol)._chopxy()
        return core.ndf(y)
    
    def f_cdf_max(self, vals, ndf):
        k = 2
        stat = (ndf - k)/k * (1-vals/vals.max())
        cdf = stats.f.cdf(stat, k, ndf)
        return cdf
        
    def f_cdf_min(self, vals, ndf):
        k = 2
        stat = (ndf - k)/k * (vals/vals.min()-1)
        cdf = stats.f.cdf(stat, k, ndf)
        return cdf
    
    @property
    def flam2pdf(self, pol=None):
        vals = self.lam2
        fast, lag = self._fast_lag_minloc(vals)
        ndf = self.ndf(fast, lag)
        return self.f_cdf_min(vals, ndf)        
    
    @property
    def flamratpdf(self, pol=None):
        vals = self.lamrat
        fast, lag = self._fast_lag_maxloc(vals)
        ndf = self.ndf(fast, lag)
        return self.f_cdf_max(vals, ndf)
    
    @property    
    def fzrhopdf(self, pol=None):
        vals = np.arctanh(self.rho)
        fast, lag = self._fast_lag_maxloc(vals)
        ndf = self.ndf(fast, lag)
        return self.f_cdf_max(vals, ndf)
        

        
    #     likelihood = 1 -f_cdf
    #     pdf = likelihood / np.sum(likelihood)
    #
    # def fresult(self, alpha=0.05):
    #     vals = self.lamrat
    #     fast, lag = self._fast_lag_maxloc(vals)
        
    def ftest_min(self, vals, ndf, alpha=0.05):
        """
        returns value (in vals) at 100(1-alpha)% confidence interval,
        by default alpha = 0.05, i.e. 95% confidence interval,
        following Silver and Chan (1991).
        """    
        k = 2 # two parameters, phi and dt.
        F = stats.f.ppf(1-alpha, k, ndf)
        val_at_alpha = vals.min() * ( 1 + (k/(ndf-k)) * F)
        return val_at_alpha
        
    def ftest_max(self, vals, ndf, alpha=0.05):
        """
        returns value (in vals) at 100(1-alpha)% confidence interval,
        by default alpha = 0.05, i.e. 95% confidence interval,
        following Silver and Chan (1991).
        """    
        k = 2 # two parameters, phi and dt.
        F = stats.f.ppf(1-alpha, k, ndf)
        val_at_alpha = vals.max() / ( 1 + (k/(ndf-k)) * F)
        return val_at_alpha        
    

    
    def ferror_min(self, vals, alpha):
        fast, lag = self._fast_lag_minloc(vals)
        ndf = self.ndf(fast, lag)
        ftestalpha = self.ftest_min(vals, ndf, alpha)
        dfast, dlag = core.contour_halfwidth(vals, ftestalpha, surftype='min')
        return fast, dfast, lag, dlag
        
    def ferror_max(self, vals, alpha):
        # CHECK THIS MATHS I JUST FILLED SOMETHING IN
         fast, lag = self._fast_lag_maxloc(vals)
         ndf = self.ndf(fast, lag)
         ftestalpha = self.ftest_max(vals, ndf, alpha)
         dfast, dlag = core.contour_halfwidth(vals, ftestalpha, surftype='max')
         return fast, dfast, lag, dlag   



    # def fast_lag_dfast_dlag(self, pdf, alpha=0.05):
        
        
    
    # def _bootstrap_kdes(self, fast, lag, n=2000, **kwargs):
    #     """Estimate distributions for pearson's r, the lam1/lam2 ratio
    #     (and, if pol is not set, the source polarisation) for the data corrected
    #     using the parameters *fast* and *lag*. The distributions are estimated using bootstrapping.
    #
    #     args (required):
    #
    #     (splitting parameters to correct data by)
    #
    #     fast -- the fast direction
    #     lag  -- the delay time
    #
    #     kwargs (optional):
    #
    #     n    -- the number of bootstrap iterations (default=2000)
    #     pol -- if pol is specified uses transverse minimisation method
    #            and makes no attempt to calculate spol.
    #     """
    #
    #     if self._rcvcorr is not None:
    #         raise NotImplementedError('Not yet implemented.')
    #
    #     if self._srccorr is not None:
    #         raise NotImplementedError('Not yet implemented.')
    #
    #     # prepare data by backing off the splitting parameters
    #     # uses wraparound in time window to replicate behaviour in frequency domain.
    #     bsindata = self.data._wrap_unsplit(fast, lag, taper=self._taper, **kwargs)
    #     x, y = bsindata.x, bsindata.y
    #
    #     # calculate bootstrap covariance matrices
    #     bscov = np.empty((n, 2, 2))
    #     for ii in range(n):
    #         bsx, bsy = core.bootstrap_resamp(x, y)
    #         bscov[ii] = core.cov2d(bsx, bsy)
    #
    #     # calculate rho kde
    #     stdx = np.sqrt(bscov[:, 0, 0])
    #     stdy = np.sqrt(bscov[:, 1, 1])
    #     rho = bscov[:, 0, 1] / (stdx * stdy)
    #     r_kde = core.kde(np.abs(rho))
    #
    #
    #     if self._pol is None:
    #         # Use eigenvalue method
    #         # rotate to zero
    #         rot = core._rot(-fast)
    #         bscov = np.matmul(rot, np.matmul(bscov, rot.T))
    #         # use eigenvector method.
    #         evals, evecs = np.linalg.eigh(bscov)
    #         rat = evals[:,1]/evals[:,0]
    #         spol = (np.rad2deg(np.arctan2(evecs[:,1,1], evecs[:,0,1]))+3690)%180-90
    #         rat_kde = core.kde(rat)
    #         spol_kde = core.kde(spol)
    #         # all done
    #         return r_kde, rat_kde, spol_kde
    #     else:
    #         # Use transverse energy method
    #         # rotate to pol
    #         rot = core._rot(pol-fast)
    #         bscov = np.matmul(rot, np.matmul(bscov, rot.T))
    #         rat = bscov[:,0,0]/bscov[:,1,1]
    #         rat_kde = core.kde(rat)
    #         return r_kde, rat_kde

        
        
    
    # def silver_chan(self, **kwargs):
    #
    #     if self._pol is None:
    #         vals = self.lam2
    #     else:
    #         vals = self.transenergy
    #
    #     return Method(self, vals=vals, ftest=True)


        
        # if self._bootstrap:
        #     bscov = self._bootcov(m.fast, m.lag, **kwargs)
        #     bsvals = core.bscov2eigrat(bscov, **kwargs)
        #     norm = core.norm(bsvals)
        #
        # m.likelihood = norm.pdf(vals)
        # m.pdf = m.likelihood / np.sum(m.likelihood)
            
            # _, m.kde, m.spol_kde = self._bootstrap_kdes(m.fast, m.lag, **kwargs)
            # m.likelihood = m.kde.pdf(vals.flatten()).reshape((vals.shape))
            # m.loglike = m.kde.logpdf(vals.flatten()).reshape((vals.shape))
            # m.pdf = m.likelihood / np.sum(m.likelihood)
            
        #
        # # error estimation
        # n = m.ndf
        # k = 2
        # stat = (n - k)/k * (lam2/lam2.min() - 1)
        # m.f_cdf = stats.f.cdf(stat, k, n)
        #
        # # calculate width of 1 sigma errors
        # critval = 0.683
        # m.dfast, m.dlag = self._contour_halfwidth(m.f_cdf, critval, surftype='min')
        
            
        # return m
        # sc = {}
        # sc['lam1'] = lam1
        # sc['lam2'] = lam2
        # sc['maxidx'] = core.max_idx(lam1/lam2)
        # dd, ll = self._grid
        # sc['fast'] = dd[sc['maxidx']]
        # sc['lag']  = ll[sc['maxidx']]
 
        # sc['srcpol'] = self._covmap[ml]
        # sc['dfast']
        # sc['dlag']
        # self.__silver_chan = sc
        
    # def correlation(self, **kwargs):
    #
    #     vals = np.abs(core.covmap2rho(self._covmap))
    #     fvals = np.arctanh(vals)
    #     return Method(self, vals=fvals)
        
        #
        # # if self._bootstrap:
        # #     m.kde, _, _ = self._bootstrap_kdes(m.fast, m.lag, **kwargs)
        # #     m.likelihood = m.kde.pdf(vals.flatten()).reshape((vals.shape))
        # #     m.loglike = m.kde.logpdf(vals.flatten()).reshape((vals.shape))
        # #     m.pdf = m.likelihood / np.sum(m.likelihood)
        #
        # if self._bootstrap:
        #     bscov = self._bootcov(m.fast, m.lag, **kwargs)
        #     bsvals = core.bscov2rho(bscov, **kwargs)
        #
        # fbsvals = np.arctanh(bsvals)
        # norm = core.norm(fbsvals)
        # m.likelihood = norm.pdf(fvals)
        # m.pdf = m.likelihood / np.sum(m.likelihood)
        #
        # # xc = {}
        # # xc['rho'] = np.abs(core.covmap2rho(self._covmap))
        # # xc['maxidx'] = core.max_idx(xc['rho'])
        # # dd, ll = self._grid
        # # xc['fast'] = dd[xc['maxidx']]
        # # xc['lag']  = ll[xc['maxidx']]
        # return m


                    
    # METHODS 
    #--------    

    
    # def srcpol(self):
    #     # recover source polarisation
    #     if self.data._pol is not None:
    #         return self.data._pol
    #     else:
    #         return self.data_corr()._estimate_pol()
    #
    # def snr(self):
    #     """Restivo and Helffrich (1999) signal to noise ratio"""
    #     d = self.srcpoldata_corr()._chop()
    #     return core.snrRH(d.x, d.y)
    #

            
    # F-test utilities
    
    # def ndf(self):
    #     """Number of degrees of freedom."""
    #     x, y = self.srcpoldata_corr()._chopdata()
    #     return core.ndf(y)
    
    # def get_errors(self, surftype=None):
    #     """
    #     Return dfast and dtlag.
    #
    #     These errors correspond to one sigma in the parameter estimate.
    #
    #     Calculated by taking a quarter of the width of 95% confidence region (found using F-test).
    #     """
    #
    #     # search interval steps
    #     lag_step = self.lags[1] - self.lags[0]
    #     fast_step = self.degs[1] - self.degs[0]
    #
    #     # Find nodes where we fall within the 95% confidence region
    #
    #     if surftype == 'max':
    #         confbool = self.errsurf >= self.conf95level
    #     elif surftype == 'min':
    #         confbool = self.errsurf <= self.conf95level
    #     else:
    #         raise ValueError('surftype must be \'min\' or \'max\'')
    #
    #     # tlag error
    #     lagbool = confbool.any(axis=1)
    #     # last true value - first true value
    #     truth = np.where(lagbool)[0]
    #     fdlag = (truth[-1] - truth[0] + 1) * lag_step * 0.25
    #
    #     # fast error
    #     fastbool = confbool.any(axis=0)
    #     # trickier to handle due to cyclicity of angles
    #     # search for the longest continuous line of False values
    #     cyclic = np.hstack((fastbool, fastbool))
    #     lengthFalse = np.diff(np.where(cyclic)).max() - 1
    #     # shortest line that contains ALL true values is then:
    #     lengthTrue = fastbool.size - lengthFalse
    #     fdfast = lengthTrue * fast_step * 0.25
    #
    #     # return
    #     return fdfast, fdlag
        
    # def _contour_halfwidth(self, surf, critval, surftype=None):
    #     """
    #     Return half width of contour (dfast, dlag) for surface *surf* at value *critval*.
    #
    #     Some common critical values:
    #     1 sigma = 0.683
    #     2 sigma = 0.954
    #     3 sigma = 0.997
    #     """
    #
    #     # search interval steps
    #     lag_step = self._lags[1] - self._lags[0]
    #     fast_step = self._degs[1] - self._degs[0]
    #
    #     # Find nodes where we fall within the 95% confidence region
    #
    #     if surftype == 'max':
    #         confbool = surf >= critval
    #     elif surftype == 'min':
    #         confbool = surf <= critval
    #     else:
    #         raise ValueError('surftype must be \'min\' or \'max\'')
    #
    #     # tlag error
    #     lagbool = confbool.any(axis=1)
    #     # last true value - first true value
    #     truth = np.where(lagbool)[0]
    #     fdlag = (truth[-1] - truth[0] + 1) * lag_step * 0.5
    #
    #     # fast error
    #     fastbool = confbool.any(axis=0)
    #     # trickier to handle due to cyclicity of angles
    #     # search for the longest continuous line of False values
    #     cyclic = np.hstack((fastbool, fastbool))
    #     lengthFalse = np.diff(np.where(cyclic)).max() - 1
    #     # shortest line that contains ALL true values is then:
    #     lengthTrue = fastbool.size - lengthFalse
    #     fdfast = lengthTrue * fast_step * 0.5
    #
    #     # return
    #     return fdfast, fdlag
        
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
                return data_corr._chopdata().rotateto(self.fast)
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
        bs = self.data_corr()
        origang = bs.cmpangs()[0]
        # replace noise sequence
        bs.rotateto(self.srcpol())
        bs.y = core.resample_noise(bs.y)
        bs.rotateto(origang)
        # reapply splitting
        # src side correction
        if self._srccorr is not None: bs = bs._split(*self._srccorr)
        # target layer correction
        bs = bs._split(self.fast, self.lag)
        # rcv side correction
        if self._rcvcorr is not None: bs = bs._split(*self._rcvcorr)
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
        Save Measurement for future referral
        """
        io.save(self, filename)
        
    def load(self, filename):
        return io.load(self, filename)

    def copy(self):
        return io.copy(self)
        
    # spit out the answer
        
    # def report(self, **kwargs):
    #     """Prints fast, lag, dfast, dlag to screen/stdout."""
    #     print('fast'.rjust(10), 'dfast'.rjust(9), 'lag'.rjust(9), 'dlag'.rjust(9))
    #     print('{0:10.2f}{1:10.2f}{2:10.2f}{3:10.2f}'.format(self.fast, self.dfast, self.lag, self.dlag))
            
    
    # Plotting
    
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
           kwargs['vals'] = self.lamrat
           kwargs['title'] = r'$\lambda_1 / \lambda_2$'
        
        self._plot(**kwargs)
    
    def _plot(self, pol=None, **kwargs):
        
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


        orig = self.srcpoldata().chop()
        fast, _, lag, _ = self.silver_chan()
        corr = self.srcpoldata_corr(fast, lag).chop()
                
        # get axis scaling
        lim = np.abs(corr.data).max() * 1.1
        ylim = [-lim, lim]
        
        # long window data
        self.data._ptr(ax0, ylim=ylim, **kwargs)

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
        laggrid, deggrid = self._grid
        cax = ax.contourf(laggrid, deggrid, kwargs['vals'], 26, cmap=kwargs['cmap'])
        cbar = plt.colorbar(cax)
        ax.set_ylabel(r'Fast Direction ($^\circ$)')
        # ax.set_xlabel('Delay Time (' + self.units + ')')
        
        # confidence region
        # if 'conf95' in kwargs and kwargs['conf95'] == True:
        #     ax.contour(self.lagmap, self.degmap, self.errsurf, levels=[self.conf95level],
        #             colors='r', alpha=.5, linewidths=3)
            
        # marker
        fast, dfast, lag, dlag = self.silver_chan()
        if 'marker' in kwargs and kwargs['marker'] == True:
            ax.errorbar(lag, fast, xerr=dlag, yerr=dfast)

        ax.set_xlim([laggrid[0,0], laggrid[-1,0]])
        ax.set_ylim([deggrid[0,0], deggrid[0,-1]])
    
        # optional title
        if 'title' in kwargs:
            ax.set_title(kwargs['title']) 
            
        # add info in text box
        if 'info' in kwargs and kwargs['info'] == True:
            textstr = '$\phi=%.1f\pm%.1f$\n$\delta t=%.2f\pm%.2f$'%\
                        (fast, dfast, lag, dlag)
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
        

# class Method(Meas):
#
#     def __init__(self, meas, **kwargs):
#
#         self.meas = meas
#
#
#         # Default Settings (These can be overidden using keywords)
#         settings = {}
#         settings['vals'] = self.meas.lamrat
#         settings['ftest'] = False
#         settings['bootstrap'] = True
#         settings['alpha'] = 0.05
#
#
#
#
#     # def maxloc(self, vals):
#     #     return core.max_idx(vals)
#     #
#     # def minloc(self, vals):
#     #     return core.min_idx(vals)
#
#     def _fast_lag_maxloc(self, vals):
#         idx = core.max_idx(vals)
#         ll, dd = self.meas_grid
#         return dd[idx], ll[idx]
#
#     def _fast_lag_minloc(self, vals):
#         idx = core.min_idx(vals)
#         ll, dd = self.meas_grid
#         return dd[idx], ll[idx]
        

    # @property
    # def fast(self):
    #     ll, dd = self.meas._grid
    #     return dd[self.maxloc]
    #
    # @property
    # def lag(self):
    #     ll, dd = self.meas._grid
    #     return ll[self.maxloc]
        
    # error estimation
    
    # @property
    # def ndf(self):
    #     x, y = self.srcpoldata_corr._chopxy()
    #     return core.ndf(y)
    #
    # # @property
    # # def dfast(self):
    # #     return 0
    # #
    # # @property
    # # def dlag(self):
    # #     return 0
    #
    # # data views
    #
    # @property
    # def data(self):
    #     return self.meas.data
    #
    # @property
    # def srcpol(self):
    #     return self.data_corr._estimate_pol()
    #
    # @property
    # def data_corr(self):
    #     data_corr = self.data
    #     # rcv side correction
    #     if self.py._rcvcorr is not None:
    #         data_corr = data_corr.unsplit(*self.py._rcvcorr)
    #     # target layer correction
    #     data_corr = data_corr.unsplit(self.fast, self.lag)
    #     # src side correction
    #     if self.py._srccorr is not None:
    #         data_corr = data_corr.unsplit(*self.py._srccorr)
    #     return data_corr
    #
    #
    #

        
    # def bscov(self, fast, lag, **kwargs):
    #     return self.meas._bootcov(fast, lag)
    #
    # def bslamrat(self, **kwargs):
    #     vals = self.meas.lamrat
    #     fast, lag = self.meas._fast_lag_maxloc(vals)
    #     return core.bscov2eigrat(self.bscov(fast, lag, **kwargs))
    #
    # def bslam2(self, **kwargs):
    #     vals = self.meas.lam2
    #     fast, lag = self.meas._fast_lag_minloc(vals)
    #     return core.bscov2lam2(self.bscov(fast, lag, **kwargs))
    #
    # def bszrho(self, **kwargs):
    #     vals = self.meas.zrho
    #     fast, slow = self.meas._fast_lag_maxloc(vals)
    #     return core.bscov2zrho(self.bscov(fast, lag, **kwargs))
    #
    #
    # def bsnormpdf(self, bsvals, vals):
    #     norm = core.norm(bsvals)
    #     likelihood = norm.pdf(vals)
    #     pdf = likelihood / np.sum(likelihood)
    #     return pdf
    #
    # def bskdepdf(self, bsvals, vals):
    #     kde = core.kde(bsvals)
    #     likelihood = kde.pdf(vals.flatten()).reshape((vals.shape))
    #     pdf = likelihood / np.sum(likelihood)
    #     return pdf
        
    # def fsurf(self):
    #     vals = self.lam2
    #     fast, lag = self.meas._fast_lag_minloc(vals)
        # ndf =       
    #
    # # plotting
    # def psurf(self, **kwargs):
    #     fig, ax = plt.subplots(1)
    #     self._psurf(ax, **kwargs)
    #     plt.show()
    #
    #
    # def _psurf(self, ax, **kwargs):
    #     """
    #     Plot an error surface.
    #
    #     **kwargs
    #     - cmap = 'magma'
    #     - vals = (M.lam1-M.lam2) / M.lam2Split
    #     """
    #
    #     if 'cmap' not in kwargs:
    #         kwargs['cmap'] = 'magma'
    #
    #     ll, dd = self.py._grid
    #
    #
    #     # error surface
    #     cax = ax.contourf(ll, dd, self.vals, 26, cmap=kwargs['cmap'])
    #     cbar = plt.colorbar(cax)
    #     ax.set_ylabel(r'Fast Direction ($^\circ$)')
    #     ax.set_xlabel('Delay Time (' + self.py.data.units + ')')
    #
    #     # confidence region
    #     # if 'conf95' in kwargs and kwargs['conf95'] == True:
    #     #     ax.contour(self.lagmap, self.degmap, self.errsurf, levels=[self.conf95level],
    #     #             colors='r', alpha=.5, linewidths=3)
    #
    #     # marker
    #     if 'marker' in kwargs and kwargs['marker'] == True:
    #         ax.errorbar(self.lag, self.fast, xerr=self.dlag, yerr=self.dfast)
    #
    #     ax.set_xlim([ll[0,0], ll[-1,0]])
    #     ax.set_ylim([dd[0,0], dd[0,-1]])
    #
    #     # optional title
    #     if 'title' in kwargs:
    #         ax.set_title(kwargs['title'])
    #
    #     # add info in text box
    #     if 'info' in kwargs and kwargs['info'] == True:
    #         textstr = '$\phi=%.1f\pm%.1f$\n$\delta t=%.2f\pm%.2f$'%\
    #                     (self.fast, self.dfast, self.lag, self.dlag)
    #         # place a text box in upper left in axes coords
    #         props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    #         ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=12,
    #                 verticalalignment='top', bbox=props)
    #
    #     # plot particle motions on surface
    #     if 'ppm' in kwargs and kwargs['ppm'] is True:
    #         sublags = self.lags[0:-1:int(self.lags.size/6)]
    #         subdegs = self.degs[0:-1:int(self.degs.size/6)]
    #         sublags = sublags + (self.lags[-1]-sublags[-1]) / 2
    #         subdegs = subdegs + (self.degs[-1]-subdegs[-1]) / 2
    #         x, y = self.SplitWave_corr()._chopdata()
    #         lagtot = self.lags[-1] - self.lags[0]
    #         degtot = self.degs[-1] - self.degs[0]
    #         boost = 10 * lagtot / np.max((x**2 + y**2)**.5)
    #         for fast in subdegs:
    #             for lag in sublags:
    #                 x, y = self.unsplit(fast, lag)._chopdata()
    #                 ax.plot(lag + y*boost/degtot, fast + x*boost/lagtot, color='w',alpha=0.5)
    #
    #     return ax
    #
    # def _plot(self, **kwargs):
    #
    #     # if 'vals' not in kwargs:
    #     #     raise Exception('vals must be specified')
    #
    #     # if kwargs['vals'] == 'pdf':
    #     #     kwargs['vals'] = self.estimate_pdf()
    #     #     kwargs['title'] = r'Probability Density'
    #
    #     # setup figure and subplots
    #     fig = plt.figure(figsize=(12,6))
    #
    #     gs = gridspec.GridSpec(3, 3,
    #                        width_ratios=[2,1,3]
    #                        )
    #     ax0 = plt.subplot(gs[0,0:2])
    #     ax1 = plt.subplot(gs[1,0])
    #     ax2 = plt.subplot(gs[1,1])
    #     ax3 = plt.subplot(gs[2,0])
    #     ax4 = plt.subplot(gs[2,1])
    #     ax5 = plt.subplot(gs[:,2])
    #
    #     orig = self.srcpoldata.chop()
    #     corr = self.srcpoldata_corr.chop()
    #
    #     # get axis scaling
    #     lim = np.abs(corr.data).max() * 1.1
    #     ylim = [-lim, lim]
    #
    #     # long window data
    #     self.data._ptr(ax0, ylim=ylim, **kwargs)
    #
    #     # original
    #     orig._ptr(ax1, ylim=ylim, **kwargs)
    #     orig._ppm(ax2, lims=ylim, **kwargs)
    #
    #     # corrected
    #     corr._ptr(ax3, ylim=ylim, **kwargs)
    #     corr._ppm(ax4, lims=ylim, **kwargs)
    #
    #     # add marker and info box by default
    #     if 'marker' not in kwargs: kwargs['marker'] = True
    #     if 'info' not in kwargs: kwargs['info'] = True
    #     if 'conf95' not in kwargs: kwargs['conf95'] = True
    #     self._psurf(ax5,**kwargs)
    #
    #     # title
    #     if 'name' in kwargs:
    #         plt.suptitle(kwargs['name'])
    #
    #     # neaten
    #     plt.tight_layout()
    #
    #     # save or show
    #     if 'file' in kwargs:
    #         plt.savefig(kwargs['file'])
    #     else:
    #         plt.show()
            
# class Pdf(Method):
#
#     def __init__(self, meth, **kwargs):
#
#         self.meth = meth
#         self.meas = meth.meas
#
#         # Default Settings (These can be overidden using keywords)
#         settings = {}
#         settings['vals'] = self.meth.meas.lamrat
#         settings['ftest'] = False
#         settings['bootstrap'] = True
#         settings['alpha'] = 0.05
#
#
#         # supported vals
#         # if self.vals not in [ self.lamrat, self.lam1, self.lam2, self.rad1, self.rad2, self.radrat]
#         # workout what function to use
#
#         def bscov(self, fast, lag, **kwargs):
#             return self.meas._bootcov(fast, lag)
#
#         def bslamrat(self, **kwargs):
#             vals = self.meas.lamrat
#             fast, lag = self.meas._fast_lag_maxloc(vals)
#             return core.bscov2eigrat(self.bscov(fast, lag, **kwargs))
#
#         def bslam2(self, **kwargs):
#             vals = self.meas.lam2
#             fast, lag = self.meas._fast_lag_minloc(vals)
#             return core.bscov2lam2(self.bscov(fast, lag, **kwargs))
#
#         def bszrho(self, **kwargs):
#             vals = self.meas.zrho
#             fast, slow = self.meas._fast_lag_maxloc(vals)
#             return core.bscov2zrho(self.bscov(fast, lag, **kwargs))
#
#
#         def bsnormpdf(self, bsvals, vals):
#             norm = core.norm(bsvals)
#             likelihood = norm.pdf(vals)
#             pdf = likelihood / np.sum(likelihood)
#             return pdf
#
#         def bskdepdf(self, bsvals, vals):
#             kde = core.kde(bsvals)
#             likelihood = kde.pdf(vals.flatten()).reshape((vals.shape))
#             pdf = likelihood / np.sum(likelihood)
#             return pdf
#
#         def bs_pdf(self, bsvals, vals, kde=False, maxnode=True, **kwargs):
#
#             # vals = self.meas.lamrat
#             # if maxnode:
#             #     fast, slow = self.meas._fast_lag_maxloc(vals)
#             # else:
#             #     fast, slow = self.
#             # bsvals = core.bscov2eigrat(self.bscov(fast, slow, **kwargs))
#
#             if kde:
#                 # fit KDE
#                 return self.bskdepdf(bsvals, vals)
#             else:
#                 # fit to normal distribution
#                 return self.bsnormpdf(bsvals, vals)
#
#
#         def bs_lamrat(self):
#             vals = self.meas.lamrat
#             bsvals = core.bscov2eigrat()
#
#         @property
#         def ftest(self, **kwargs):
#             return
    