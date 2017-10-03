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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class EigenM:
    
    """
    Silver and Chan (1991) eigenvalue method measurement.
    
    >>> EigenM(**kwargs) returns a synthetic measurement
    
    >>> EigenM(Pair,**kwargs) returns a measurement on data contained in Pair
    
    kwargs:
    rcvcorr = (fast,tlag) | tuple
    srccorr = (fast,tlag) | tuple
    
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
        
        # other stuff
        self.units = self.data.units
        
        # process kwargs
              
        if ('tlags' in kwargs):
            # round to nearest 2
            lags = 2 * np.rint( 0.5 * kwargs['tlags'] / self.delta )
            kwargs['lags'] = np.unique(lags).astype(int)
                      
        if ('rcvcorr' in kwargs):
            # convert time shift to nsamples -- must be even
            degree = kwargs['rcvcorr'][0]
            nsamps = int(2 * np.rint( 0.5 * kwargs['rcvcorr'][1] / self.delta ))
            kwargs['rcvcorr'] = (degree, nsamps)
            self.rcvcorr = ( degree, nsamps * self.delta)
                   
        if ('srccorr' in kwargs):
            # convert time shift to nsamples -- must be even
            degree = kwargs['srccorr'][0]
            nsamps = int(2 * np.rint( 0.5 * kwargs['srccorr'][1] / self.delta ))
            kwargs['srccorr'] = (degree, nsamps)
            self.srccorr = (degree,nsamps * self.delta)

            
        # ensure trace1 at zero angle
        self.data.rotateto(0)        
        
        # grid search splitting
        self.degs, self.lags, self.lam1, self.lam2, self.window = eigval.grideigval(self.data.x,self.data.y,**kwargs)
            
        self.tlags = self.lags * self.delta
                
        # get some measurement attributes
        # uses ratio lam1/lam2 to find optimal fast and lag parameters
        maxloc = core.max_idx(self.lam1/self.lam2)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        self.tlag = self.lag * self.delta

        
        # generate "squashed" profiles
        self.fastprofile = np.sum(self.lam1/self.lam2, axis=0)
        self.lagprofile = np.sum(self.lam1/self.lam2, axis=1)
        # generate redefined "NI" value
        self.ni = ni(self)
        
        # correct the data
        x,y = self.data.x, self.data.y
        # rcv side      
        if 'rcvcorr' in kwargs:
            x,y = core.unsplit(x,y,*kwargs['rcvcorr'])
        # target layer
        x,y = core.unsplit(x,y,self.fast,self.lag)
        # src side
        if 'srccorr' in kwargs:
            x,y = core.unsplit(x,y,*kwargs['srccorr'])
        
        self.data_corr = Pair(x,y,**kwargs)
        
        self.srcpol = core.pca(self.data_corr.x,self.data_corr.y)
        
        # signal to noise calculation
        ### if total energy = signal + noise = lam1 + lam2
        ### lam1 = signal + 1/2 noise
        ### lam2 = 1/2 noise
        ### then signal = lam1 - lam2
        ###       noise = 2 * lam2        
        self.snr = np.max((self.lam1-self.lam2)/(2*self.lam2))
        
        # error estimations
        self.dfast, self.dtlag = self.f_errors()
    
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
        tlag_step = self.tlags[1,0]-self.tlags[0,0]
        fast_step = self.degs[0,1]-self.degs[0,0]
        
        # Find nodes where we fall within the 95% confidence region
        confbool = self.lam2 <= self.lam2_95() 
        
        # tlag error
        tlagbool = confbool.any(axis=1)
        # last true value - first true value
        truth = np.where(tlagbool)[0]
        dtlag = (truth[-1] - truth[0] + 1) * tlag_step * 0.25   
          
        # fast error
        fastbool = confbool.any(axis=0)
        # trickier to handle due to cyclicity of angles
        # search for the longest continuous line of False values
        cyclic = np.hstack((fastbool,fastbool))
        lengthFalse = np.diff(np.where(cyclic)).max()
        # shortest line that contains ALL true values is then:
        lengthTrue = fastbool.size - lengthFalse
        dfast = lengthTrue * fast_step * 0.25 
               
        # return
        return dfast, dtlag
        
    # Output
    
    def report(self):
        """
        Print out a summary of the result.
        """
        
        # pformat = 

    # def save():
    #     """
    #     Save Measurement for future referral
    #     """     
        
    
    # Plotting

    def plot(self,**kwargs):
        
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
        
        # original
        plot.trace(d1.x,d1.y,time=d1.t(),ax=ax0)
        plot.particle(d1.x,d1.y,ax=ax1)
        
        # corrected
        plot.trace(d2.x,d2.y,time=d2.t(),ax=ax2)
        plot.particle(d2.x,d2.y,ax=ax3)
        
        # error surface
        plot.surf(self,ax=ax4,**kwargs)
        
        # neaten and show
        plt.tight_layout()
        plt.show()


def ni(M):
    """
    measure of self-similarity in measurements at 90 degree shift in fast direction
    """
    halfway = int(M.degs.shape[1]/2)
    diff = M.fastprofile - np.roll(M.fastprofile,halfway)
    mult = M.fastprofile * np.roll(M.fastprofile,halfway)
    sumdiffsq = np.sum(diff**2)
    summult = np.sum(mult)
    return sumdiffsq/summult