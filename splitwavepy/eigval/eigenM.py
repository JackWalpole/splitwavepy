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

import numpy as np
import matplotlib.pyplot as plt

# for contouring
# from skimage import measure

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
    
    def srcpoldata(self):
        return Pair(*core.rotate(self.data.x,self.data.y,self.srcpol))
        
    def srcpoldata_corr(self):
        data_corr = self.data_corr
        return Pair(*core.rotate(data_corr.x,data_corr.y,self.srcpol))
    
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
        
    # def dfast(self):
    #     """
    #     Standard error in fast direction using F test.
    #
    #     Quarter width of 95% confidence contour along fast axis.
    #     """
    #
    #
    # def dtlag(self):
    #     """
    #     Standard error in delay time using F test.
    #
    #     Quarter width of 95% confidence contour along lag axis.
    #     """
    
    # Output
    
    def report(self):
        """
        Print out a summary of the result.
        """
        
    
    # Plotting
    
    def plotsurf(self,vals=None,cmap='magma',lam2_95=False,polar=False):
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
                plt.contour(rads,lags,lam2,levels=[self.lam2_95()])
        else:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            # fig = plt.figure(figsize=(6,6)) 
            
            ax = plt.subplot(111)
            
            # error surface
            # v = np.linspace(0, 50, 26, endpoint=True)
            # cax = ax.contourf(self.tlags,self.degs,vals,v,cmap=cmap,extend='max')
            cax = ax.contourf(self.tlags,self.degs,vals,26,cmap=cmap,extend='max')
            ax.set_yticks(np.linspace(-90,90,7,endpoint=True))
            # cbar = plt.colorbar(cax,ticks=v[::5])
            cbar = plt.colorbar(cax)
            marker = plt.plot(self.tlag,self.fast,'k+',markersize=10.)
                   
            if lam2_95 is True:
                plt.contour(self.tlags,self.degs,self.lam2,levels=[self.lam2_95()])
            
            
            # create new axes on the right and on the top of the current axes.
            divider = make_axes_locatable(ax)

            axFast = divider.append_axes("left", size=.5, pad=0.4, sharey=ax)
            axFast.plot(self.fastprofile,self.degs[0,:])
            axFast.axes.get_xaxis().set_visible(False)
            axFast.axes.get_yaxis().set_visible(False)
            axFast.set_title(r'Fast Direction (degrees)')
            axFast.invert_xaxis()
            # axFast.fill_betweenx(color='b',alpha=0.2)
            
            # ax.set_xlabel(r'Delay Time (s)')
            # ax.set_ylabel(r'Fast Direction (degrees)')
            
            axLag  = divider.append_axes("bottom", size=.5, pad=0.4, sharex=ax)
            axLag.plot(self.tlags[:,0],self.lagprofile)
            axLag.axes.get_xaxis().set_visible(False)
            axLag.axes.get_yaxis().set_visible(False)
            axLag.invert_yaxis()
            
            ax.set_xlim(self.tlags[0][0],self.tlags[-1][-1])
            ax.set_ylim(self.degs[0][0],self.degs[-1][-1])
        
        plt.show()

    # def save():
    #     """
    #     Save Measurement for future referral
    #     """
    
    def plot(M):
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(12,6)) 
        gs = gridspec.GridSpec(2, 3,
                           width_ratios=[1,1,2]
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
        # ax1 -- trace orig
        ax1.plot(d1.t(),d1.x)
        ax1.plot(d1.t(),d1.y)
        ax1.axes.get_yaxis().set_visible(False)
        # ax2 -- hodo orig
        lim = np.abs(d1.xy()).max() * 1.1
        ax2.axis('equal')
        ax2.plot(d1.y,d1.x)
        ax2.set_xlim([-lim,lim])
        ax2.set_ylim([-lim,lim])
        ax2.axes.get_xaxis().set_visible(False)
        ax2.axes.get_yaxis().set_visible(False)
        # ax3 -- trace new
        ax3.plot(d2.t(),d2.x)
        ax3.plot(d2.t(),d2.y)
        ax3.axes.get_yaxis().set_visible(False)
        # ax4 -- hodo new
        lim = np.abs(d2.xy()).max() * 1.1
        ax4.axis('equal')
        ax4.plot(d2.y,d2.x)
        ax4.set_xlim([-lim,lim])
        ax4.set_ylim([-lim,lim])
        ax4.axes.get_xaxis().set_visible(False)
        ax4.axes.get_yaxis().set_visible(False)
        # ax5 -- error surface
        # v = np.linspace(0, 50, 26, endpoint=True)
        cax = ax5.contourf(M.tlags,M.degs,vals,26,cmap='magma',extend='max')
        ax5.set_xlabel(r'Delay Time (s)')
        ax5.set_ylabel(r'Fast Direction (degrees)')
        cbar = plt.colorbar(cax)
        # cbar = plt.colorbar(cax,ticks=v[::5])
        
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