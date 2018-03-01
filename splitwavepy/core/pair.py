# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, geom
from .data import Data, Window, WindowPicker
# from .measure import Measure
from .eigenM import EigenM
from .xcorrM import XcorrM
from .transM import TransM

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
    def measureEigenM(self, **kwargs):        
        self.EigenM = EigenM(self.data, **kwargs)
        
    def measureXcorrM(self, **kwargs):
        self.XcorrM = XcorrM(self.data, **kwargs)
        
    def measureTransM(self, **kwargs):
        self.TransM = TransM(self.data, **kwargs)
        
    # Other 
    def splitting_intensity(self, **kwargs):
        """
        Calculate the splitting intensity as defined by Chevrot (2000).
        """
        
        if 'pol' not in kwargs:
            raise Exception('pol must be specified')
            
        copy = self.data.copy()
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

    # def data_corr(self, fast, lag, **kwargs):
    #     # copy data
    #     data_corr = self.copy()
    #     # rcv side correction
    #     if kwargs['rcvcorr'] is not None:
    #         data_corr.unsplit(*kwargs['rcvcorr'])
    #     # target layer correction
    #     data_corr.unsplit(fast, lag)
    #     # src side correction
    #     if kwargs['srccorr'] is not None:
    #         data_corr.unsplit(*kwargs['srccorr'])
    #     return data_corr

      
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

