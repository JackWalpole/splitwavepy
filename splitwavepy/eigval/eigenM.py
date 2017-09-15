"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
from ..core import pair
from ..core import window
from . import eigval

import numpy as np
import matplotlib.pyplot as plt

class EigenM:
    
    """
    Silver and Chan (1991) eigenvalue method measurement.
    """
    
    def __init__(self,*args,tlags=None,degs=None,window=None,rcvcorr=None,srccorr=None,**kwargs):
        """
        Populates an EigenM instance.
        """
        
        # process input
        if len(args) == 1 and isinstance(args[0],pair.Pair):
            self.data = args[0]
        else:
            self.data = pair.Pair(*args,**kwargs)
            
        # ensure trace1 at zero angle
        self.data.rotateto(0)
        
        # convert times to nsamples
        self.delta = self.data.delta
        
        if tlags is not None:
            lags = tlags / self.delta
        else:
            lags = None
            
        # if window is not None and isinstance(window,window.Window):
        #     self.window = window
            
        if rcvcorr is not None:
            # convert time shift to nsamples -- must be even
            nsamps = int(rcvcorr[1]/self.delta)
            nsamps = nsamps if nsamps%2==0 else nsamps + 1
            rcv = (rcvcorr[0],nsamps)
        else:
            rcv = None
        
        if srccorr is not None:
            nsamps = int(srccorr[1]/self.delta)
            nsamps = nsamps if nsamps%2==0 else nsamps + 1
            src = (srccorr[0],nsamps)
        else:
            src = None
        
        # grid search splitting
        self.degs, self.lags, self.lam1, self.lam2, self.window = eigval.grideigval(
                                                                        self.data.data,lags=lags,degs=degs,
                                                                        window=window,rcvcorr=rcv,srccorr=src)
        self.tlags = self.lags * self.delta
        
        self.rcvcorr = rcvcorr
        self.srccorr = srccorr
        
        # get some measurement attributes
        # uses ratio lam1/lam2 to find optimal fast and lag parameters
        maxloc = core.max_idx(self.lam1/self.lam2)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        # generate "squashed" profiles
        self.fastprofile = np.sum(self.lam1/self.lam2, axis=0)
        self.lagprofile = np.sum(self.lam1/self.lam2, axis=1)
        # generate redefined "NI" value
        self.ni = ni(self)
        
        # get some useful stuff
        self.data_corr = core.unsplit(self.data.data,self.fast,self.lag)
        self.srcpol = core.pca(self.data_corr)
        self.srcpoldata = core.rotate(self.data.data,-self.srcpol)
        self.srcpoldata_corr = core.rotate(self.data_corr,-self.srcpol)
        
        # signal to noise ratio estimates
        # self.snr = c.snr(c.window(self.srcpoldata_corr,self.window))
        # self.snrRH = c.snrRH(c.window(self.srcpoldata_corr,self.window))
        # self.snr = np.max(self.lam1/self.lam2)
        ### if total energy = signal + noise = lam1 + lam2
        ### lam1 = signal + 1/2 noise
        ### lam2 = 1/2 noise
        ### then signal / noise = 
        self.snr = np.max((self.lam1-self.lam2)/(2*self.lam2))

        # number degrees of freedom
        self.ndf = eigval.ndf(core.chop(self.srcpoldata_corr[1,:],self.window.width,self.window.centre,self.window.tukey))
        # value of lam2 at 95% confidence contour
        self.lam2_95 = eigval.ftest(self.lam2,self.ndf,alpha=0.05)

        # convert traces to Pair class for convenience
        self.data_corr = pair.Pair(self.data_corr,delta=self.delta)
        self.srcpoldata = pair.Pair(self.srcpoldata,delta=self.delta)
        self.srcpoldata_corr = pair.Pair(self.srcpoldata_corr,delta=self.delta)
        

    def plotsurf(self,vals=None,cmap='viridis',lam2_95=True,polar=False):
        """
        plot the measurement.
        by default plots lam1/lam2 with the lambda2 95% confidence interval overlaid
        """
        
        
        
        if vals is None:
            vals = self.lam1 / self.lam2
        
        if polar is True:
            rads = np.deg2rad(np.column_stack((self.degs,self.degs+180,self.degs[:,0]+360)))
            lags = np.column_stack((self.tlags,self.tlags,self.tlags[:,0]))
            vals = np.column_stack((vals,vals,vals[:,0]))
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            ax.contourf(rads,lags,vals,50,cmap=cmap)
            ax.set_theta_direction(-1)
            ax.set_theta_offset(np.pi/2.0)
            if lam2_95 is True:
                lam2 = np.column_stack((self.lam2,self.lam2,self.lam2[:,0]))
                plt.contour(rads,lags,lam2,levels=[self.lam2_95])
        else:
            plt.contourf(self.tlags,self.degs,vals,50,cmap=cmap)        
            if lam2_95 is True:
                plt.contour(self.tlags,self.degs,self.lam2,levels=[self.lam2_95])
            

            
        
        plt.show()

    # def save():
    #     """
    #     Save Measurement for future referral
    #     """
    
    def plot(M):
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(12, 3)) 
        gs = gridspec.GridSpec(2, 3,
                           width_ratios=[4,1,2]
                           )
    
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[1,0])
        ax4 = plt.subplot(gs[1,1])
        ax5 = plt.subplot(gs[:,2])
        
        d1 = M.data.copy()
        d1.chop(M.window)
        d2 = M.data_corr.copy()
        d2.chop(M.window)
    
        vals = M.lam1 / M.lam2
        ax1.plot(d1.t(),d1.data[0])
        ax1.plot(d1.t(),d1.data[1])
        lim = abs(d1.data.max()) * 1.1
        ax2.axis('equal')
        ax2.plot(d1.data[1],d1.data[0])
        ax2.set_xlim([-lim,lim])
        ax2.set_ylim([-lim,lim])
        ax2.axes.get_xaxis().set_visible(False)
        ax2.axes.get_yaxis().set_visible(False)
        ax3.plot(d2.t(),d2.data[0])
        ax3.plot(d2.t(),d2.data[1])
        lim = abs(d2.data.max()) * 1.1
        ax4.axis('equal')
        ax4.plot(d2.data[1],d2.data[0])
        ax4.set_xlim([-lim,lim])
        ax4.set_ylim([-lim,lim])
        ax4.axes.get_xaxis().set_visible(False)
        ax4.axes.get_yaxis().set_visible(False)
        ax5.contourf(M.tlags,M.degs,vals,20,cmap='magma')
    
        plt.show()
# def _synthM(deg=25,lag=10):
#     P = c.Pair()
#     P.split(deg,lag)
#     return eigval.grideigval(P.data)

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