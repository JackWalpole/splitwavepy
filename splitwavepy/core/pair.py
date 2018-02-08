# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, geom
from .data import Data, Window, WindowPicker
# from .measure import Measure
from .eigenM import EigenM

import pickle
import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection


class Pair:
    """
    The Pair: work with 2-component data.
        
    Usage: Pair(**kwargs)     => create Pair of synthetic data
           Pair(x,y) => creates Pair from two traces stored in numpy arrays x and y.
    
    Keyword Arguments:
        - delta = 1. (sample interval) [default] | float
        # - t0 = 0. (start time) DEVELOPMENT
    
    Naming Keyword Arguments:
        - name = 'untitled' (should be unique identifier) | string
        - cmplabels = ['cmp1','cmp2'] | list of strings
        # - units = 's' (for labelling) | string
    
    Geometry Keyword Arguments (if in doubt don't use):
        # - geom = 'geo' (N,E) / 'ray' (SV,SH) / 'cart' (X,Y)
        # - cmpvecs = np.eye(2)
        # - rcvloc = (lat,lon,r) / (x,y,z)
        # - srcloc =  (lat,lon,r) / (x,y,z)
        # - rayloc =  rcvloc # arbirary point specified by user
        # - rayvec = [0,0,1] # wave plane normal

    Methods:
        # display
            - plot() # plot traces and partice motion
            - ppm() # plot particle motion
            - ptr() # plot waveform trace
        # splitting
            - split(fast,lag)
            - unsplit(fast,lag)
        # useful
            - t()           # time
            - data()        # both x and y in one array
            - chop()        # windowed data
            - rotateto()
        # windowing
            - set_window(start,end,Tukey=None)         
        # io
            - copy()
            - save()
        # window picker
            - plot(pick=True)
    """
    def __init__(self, *args, **kwargs):
        
        # if no args make synthetic
        if len(args) == 0:
            x, y = core.synth(**kwargs)
        # otherwise read in data
        elif len(args) == 2:
            if not (isinstance(args[0], np.ndarray) & isinstance(args[1], np.ndarray)):
                raise TypeError('expecting numpy arrays')
            x, y = args[0], args[1]
        else:
            raise Exception('Unexpected number of arguments')
        
        # Initialise Data
        self.data = Data(x, y, *args, **kwargs)

    # METHODS
      



    # Measurement
    def EigenM(self, **kwargs):
        
        self.EigenM = EigenM(self.data, **kwargs)
        # M = Measure(self.data, **kwargs)
        #
        # # MAKE MEASUREMENT
        # stuff = np.asarray(M.gridsearch(core.eigvalcov, **kwargs))
        # self.lam1, self.lam2 = stuff[:,:,1].T, stuff[:,:,0].T
        # maxloc = core.max_idx(self.lam1/self.lam2)
        
    # class EigenM(Measure):
    #     """
    #     Silver and Chan (1991) eigenvalue method measurement.
    #
    #     args:
    #     None = create synthetic
    #     Pair = Measure splitting on Pair object
    #     x, y = Measure splitting on traces x, and y.
    #
    #     kwargs:
    #
    #     name -- string = 'Untitled'
    #
    #     lags -- tuple = (maxlag,)
    #          -- tuple = (maxlag,Nlags)
    #          -- tuple = (minlag,maxlag,Nlags)
    #          -- numpy ndarray
    #
    #     degs -- int = degs
    #          -- numpy ndarray
    #
    #     rcvcorr = (fast,tlag) | tuple | Receiver Correction
    #     srccorr = (fast,tlag) | tuple | Source Correction
    #     """
    #
    #     def __init__(self, *args, **kwargs):
    #         """
    #         Populates an EigenM instance.
    #         """
            #
            # # process input
            # if len(args) == 1 and isinstance(args[0],Pair):
            #     self.data = args[0]
            # else:
            #     self.data = Pair(*args,**kwargs)
            #
            # Derive from Measure

        
        
        
    #
    # def set_pol(self,*args):
    #     if len(args) == 0:
    #         self.pol = self.get_pol()
    #     elif len(args) == 1:
    #         self.pol = float(args[0])
    #     else:
    #         raise Exception('Unexpected number of arguments')
    #     return
    
    # Utility 
    
  


    # def get_pol(self):
    #     """Return principal component orientation"""
    #     # rotate to zero
    #     rot = self.cmpvecs.T
    #     data = self.chop().data()
    #     xy = np.dot(rot,data)
    #     _,eigvecs = core.eigcov(xy)
    #     x,y = eigvecs[:,0]
    #     pol = np.rad2deg(np.arctan2(y,x))
    #     return pol
        
  
    
    # def chop(self):
    #     """
    #     Chop data to window
    #     """
    #     chop = self.copy()
    #     chop.x, chop.y = core.chop(chop.x, chop.y, window=chop.window)
    #     chop.window.offset = 0
    #     return chop

    
    def splitting_intensity(self, **kwargs):
        """
        Calculate the splitting intensity as defined by Chevrot (2000).
        """
        
        if 'pol' not in kwargs:
            raise Exception('pol must be specified')
            
        copy = self.copy()
        copy.rotateto(kwargs['pol'])
        copy.x = np.gradient(copy.x)
        rdiff, trans = copy.chopdata()
        s = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
        return s

        

        
    # def grid_eigen(self, **kwargs):
    #     """Grid search for splitting parameters using the transverse energy minimisation
    #        eigenvalue method (Silver and Chan, 1991)"""
    #     # MAKE MEASUREMENT
    #     stuff = np.asarray(self._gridsearch(core.eigvalcov, **kwargs))
    #     lam1, lam2 = stuff[:,:,1].T, stuff[:,:,0].T
    #     return lam1, lam2
    #
    # def grid_trans(self, **kwargs):
    #     """Grid search for splitting parameters using the transverse energy minimisation
    #        user-specified polarisation method (Silver and Chan, 1998)"""
    #
    #     if 'pol' not in kwargs:
    #         raise Exception('pol must be specified')
    #
    #     # MAKE MEASUREMENT
    #     stuff = np.asarray(self._gridsearch(core.transenergy, **kwargs))
    #     enrgy1, enrgy2 = stuff[:,:,1].T, stuff[:,:,0].T
    #     return enrgy1, enrgy2
    #
    # def grid_xcorr(self, **kwargs):
    #     """Grid search for splitting parameters using the cross correlation method (Ando, 1980)"""
    #     # MAKE MEASUREMENT
    #     stuff = np.asarray(self._gridsearch(core.transenergy, **kwargs))
    #     xc = stuff[:,:,0].T
    #     return xc
    #
    # def eigenM(self, **kwargs):
    #
    #     # setup dictionary to hold measurement
    #     self.eigenM = {}
    #
    #     # get degs, lags and slags
    #     self.eigenM['degs'], self.eigenM['lags'], _ = self._get_degs_lags_slags(self, **kwargs)
    #     # source and receiver corrections
    #
    #
    #     # make measurement
    #     self.eigenM['lam1'], self.eigenM['lam2'] = self.grid_eigen(self, **kwargs)
    #
    #     # get useful info
    #     maxidx = core.max_idx(lam1/lam2)
    #     fast = DEGS[maxloc]
    #     tlag  = LAGS[maxloc]
    #
    #     # estimate error
    #     core.ftest(self.lam2, self.ndf(), alpha=0.05)
    #
    #     # Populate dictionary object
    #     self.eigenM = {'lags': lags, 'degs': degs,
    #                    'rcvcorr': kwargs['rcvcorr'], 'srccorr': kwargs['srccorr'],
    #                    'lam1': lam1, 'lam2': lam2, 'maxidx': maxidx,
    #                    'fast': fast, 'tlag': tlag, 'dfast': dfast, 'dtlag': dtlag
    #                    }

    def data_corr(self, fast, lag, **kwargs):
        # copy data     
        data_corr = self.copy()
        # rcv side correction     
        if kwargs['rcvcorr'] is not None:
            data_corr.unsplit(*kwargs['rcvcorr'])    
        # target layer correction
        data_corr.unsplit(fast, lag)  
        # src side correction
        if kwargs['srccorr'] is not None:
            data_corr.unsplit(*kwargs['srccorr'])
        return data_corr
                
    # Common methods    
    
    # def _gridsearch(self, func, **kwargs):
    #
    #     """
    #     Grid search for splitting parameters applied to data using the function defined in func
    #     rcvcorr = receiver correction parameters in tuple (fast,lag)
    #     srccorr = source correction parameters in tuple (fast,lag)
    #     """
    #
    #     # get degs, lags and slags
    #     degs, _, slags = self._get_degs_lags_slags(self, **kwargs)
    #
    #     # receiver correction
    #     rcvcorr = None
    #     if ('rcvcorr' in kwargs):
    #         if not isinstance(kwargs['rcvcorr'],tuple): raise TypeError('rcvcorr must be tuple')
    #         if len(kwargs['rcvcorr']) != 2: raise Exception('rcvcorr must be length 2')
    #         # convert time shift to nsamples -- must be even
    #         deg, lag = kwargs['rcvcorr']
    #         samps = core.time2samps(lag, self.delta, 'even')
    #         rcvcorr = (deg, samps)
    #
    #     # source correction
    #     srccorr = None
    #     if ('srccorr' in kwargs):
    #         if not isinstance(kwargs['srccorr'],tuple): raise TypeError('srccorr must be tuple')
    #         if len(kwargs['srccorr']) != 2: raise Exception('srccorr must be length 2')
    #         # convert time shift to nsamples -- must be even
    #         deg, lag = kwargs['srccorr']
    #         samps = core.time2samps(lag, self.delta, 'even')
    #         srccorr = (deg, samps)
    #
    #     # avoid using "dots" in loops for performance
    #     rotate = core.rotate
    #     lag = core.lag
    #     chop = core.chop
    #     unsplit = core.unsplit
    #
    #     # ensure trace1 at zero angle
    #     copy = self.copy()
    #     copy.rotateto(0)
    #     x, y = copy.x, copy.y
    #
    #     # pre-apply receiver correction
    #     if 'rcvcorr' in kwargs:
    #         rcvphi, rcvlag = rcvcorr
    #         x, y = unsplit(x, y, rcvphi, rcvlag)
    #
    #     ######################
    #     # inner loop function
    #     ######################
    #
    #     # source correction
    #
    #     if 'srccorr' in kwargs:
    #         srcphi, srclag = srccorr
    #         def srccorr(x, y, ang):
    #             x, y = unsplit(x, y, srcphi-ang, srclag)
    #             return x, y
    #     else:
    #         def srccorr(x, y, ang):
    #             return x, y
    #
    #     # rotate to polaristation (needed for tranverse min)
    #     if 'pol' in kwargs:
    #         pol = kwargs['pol']
    #         def rotpol(x, y, ang):
    #             # rotate to pol
    #             x, y = rotate(x, y, pol-ang)
    #             return x, y
    #     else:
    #         def rotpol(x, y, ang):
    #             return x, y
    #
    #     # actual inner loop function
    #     def process(x, y, ang, shift):
    #         # remove shift
    #         x, y = lag(x, y, -shift)
    #         x, y = srccorr(x, y, ang)
    #         x, y = chop(x, y, window=self.window)
    #         x, y = rotpol(x, y, ang)
    #         return func(x, y)
    #
    #     # Do the grid search
    #     prerot = [ (rotate(x, y, ang), ang) for ang in degs ]
    #
    #     out = [ [ process(data[0], data[1], ang, shift) for shift in slags ]
    #             for (data, ang) in prerot  ]
    #
    #     return out
      
    def save(self,filename):
        """
        Save me to a file
        """       
        with open(filename, 'wb') as f:
            pickle.dump(self,f)        
            
        
    # Special
    
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

