# -*- coding: utf-8 -*-
"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
from .eigenM import SC
from .xcorrM import XC
from .measure import Measure

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


class Q(SC, XC, Measure):
    
    """
    Measure splitting using 
    
    """
    
    def __init__(self, data, **kwargs): 
        
        # Use bootstrapping
        if 'bootstrap' not in kwargs:
            kwargs['bootstrap'] = True
        
        # Initiate Measurment object with dummy func
        def dummy():
            return float('NaN')
        
        Measure.__init__(self, data, dummy, **kwargs) 
        
        # Make Silver and Chan and XCorrelation Measurements
        self.sc = data.SC(**kwargs)
        self.xc = data.XC(**kwargs)
        
        # Combine Measurements
        self.q = core.q(self.sc.fast, self.sc.lag, self.xc.fast, self.xc.lag)
        self.pdf = self._pdf(**kwargs)
        self.maxloc = core.max_idx(self.pdf)
        deggrid, laggrid = self._grid_degs_lags()
        self.fast = deggrid[self.maxloc]
        self.lag  = laggrid[self.maxloc]
        self.conf95level = self.bootstrap_conf95()
        self.errsurf = self.pdf
        self.dfast, self.dlag = self.get_errors(surftype='max')
    
    def report(self, **kwargs):
        print('SCfast'.rjust(10), 'SCdfast'.rjust(9), 'SClag'.rjust(9), 'SCdlag'.rjust(9),
              'XCfast'.rjust(9), 'XCdfast'.rjust(9), 'XClag'.rjust(9), 'XCdlag'.rjust(9), 'Q'.rjust(9))
        print('{0:10.2f}{1:10.2f}{2:10.2f}{3:10.2f}{4:10.2f}{5:10.2f}{6:10.2f}{7:10.2f}{8:10.2f}'.format(
               self.sc.fast, self.sc.dfast, self.sc.lag, self.sc.dlag, 
               self.sc.fast, self.sc.dfast, self.sc.lag, self.sc.dlag, self.q))
               
    def _pdf(self, **kwargs):
        
        sc = self.sc.pdf(**kwargs)
        xc = self.xc.pdf(**kwargs)
        return (sc + xc)/2
        
    def plot(self, **kwargs):
        # error surface
        if 'vals' not in kwargs:
           kwargs['vals'] = self.pdf
           kwargs['title'] = r'PDF'
    
        self._plot(**kwargs)
        
    # def plot(self, **kwargs):
        
        # setup figure and subplots
        # fig = plt.figure(figsize=(10,6)) 
        
        # gs = gridspec.GridSpec(5, 4, width_ratios=[3,1,3,1])
        #
        #
        # # trace data at top
        # ax0 = plt.subplot(gs[0,0:3])
        # ax1 = plt.subplot(gs[0,3])
        # self.data._ptr(ax0, **kwargs)
        # self.data._ppm(ax1, **kwargs)
        #
        # # Silver and Chan Result
        # ax01 = plt.subplot(gs[1:4,0:2])
        # ax02 = plt.subplot(gs[4,0])
        # ax03 = plt.subplot(gs[4,1])
        # self.sc._psurf(ax01, vals=self.sc.estimate_pdf())
        # corr = self.sc.data_corr().chop()
        # corr._ptr(ax02)
        # corr._ppm(ax03)
        #
        # # Cross Correlation Result
        # ax11 = plt.subplot(gs[1:4,2:4])
        # ax12 = plt.subplot(gs[4,2])
        # ax13 = plt.subplot(gs[4,3])
        # self.xc._psurf(ax11, vals=self.xc.estimate_pdf())
        # corr = self.sc.data_corr().chop()
        # corr._ptr(ax12)
        # corr._ppm(ax13)
        #
        # plt.tight_layout()
        
        # plt.show()