from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core
from . import plotting
from ..eigval.eigenM import EigenM
from . import window

import numpy as np
import matplotlib.pyplot as plt
import copy

class Pair:
    """
    The Pair is a class to store two traces in the x and y directions.
    Methods are included to facilitate analysis on this Pair of traces.
    If data is not provided on initiation will return a ricker wavelet with noise.
    Usage: Pair()     => create Pair of synthetic data
           Pair(data) => creates Pair from two traces stored as rows in numpy array data
           Pair(x,y) => creates Pair from two traces stored in numpy arrays x and y.
    Optional:
        - delta = x.  Where x = sample interval.  Default x=1.0.
        - angle = x.  Where x = angle of component in Pair.data[0]. e.g. clockwise from North (or SV "up" if in ray frame).  Default is x=0.0.
    """
    def __init__(self,*args,delta=None,angle=None,**kwargs):
        
        if delta is None:
            self.delta = 1.
        else:
            self.delta = float(delta)
            
        if angle is None:
            self.angle = 0.
        else:
            self.angle = float(angle)
        
        if len(args) == 0:
            if ('lag' in kwargs):
                # convert time shift to nsamples -- must be even
                nsamps = int(kwargs['lag']/self.delta)
                nsamps = nsamps if nsamps%2==0 else nsamps + 1
                kwargs['lag'] = nsamps                                      
            self.data = core.synth(**kwargs)            
        elif len(args) == 1:       
            self.data = args[0]       
        elif len(args) == 2:            
            self.data = np.vstack((args[0],args[1]))     
        else: 
            raise Exception('Unexpected number of arguments')
                    
        # some sanity checks
        if self.data.ndim != 2:
            raise Exception('data must be two dimensional')
        if self.data.shape[0] != 2:
            raise Exception('data must contain two traces in two rows')
        if self.data.shape[1]%2 == 0:
            raise Exception('traces must have odd number of samples')
            
    # methods
    
    # set time from start
    def t(self):
        return np.arange(self.data.shape[1]) * self.delta

    def power(self):
        return self.data[0]**2+self.data[1]**2
        
    def centre(self):
        return int(self.data.shape[1]/2)

    def plot(self,window=None):
        """
        Plot trace data and particle motion
        """
        from matplotlib import gridspec
        fig = plt.figure(figsize=(12, 3)) 
        if window is None:
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
            ax0 = plt.subplot(gs[0])
            ax0.plot(self.t(),self.data[0])
            ax0.plot(self.t(),self.data[1])
            # particle  motion
            lim = abs(self.data.max()) * 1.1
            # the polar axis:
            # ax_polar = plt.subplot(gs[1], polar=True, frameon=False)
            # ax_polar.set_rmax(lim)
            # ax_polar.patch.set_facecolor(111)
            # ax_polar.get_xaxis.set_visible(False)
            # ax_polar.grid(True)
            # the data
            ax1 = plt.subplot(gs[1])
            # ax1.patch.set_alpha(0)
            ax1.axis('equal')
            ax1.plot(self.data[1],self.data[0])
            ax1.set_xlim([-lim,lim])
            ax1.set_ylim([-lim,lim])
            ax1.axes.get_xaxis().set_visible(False)
            ax1.axes.get_yaxis().set_visible(False)
        else:
            gs = gridspec.GridSpec(1, 3, width_ratios=[3,1,1])
            ax0 = plt.subplot(gs[0])
            ax0.plot(self.t(),self.data[0])
            ax0.plot(self.t(),self.data[1])
            # the window limits
            nsamps = self.t().size
            wbeg = window.start(nsamps)*self.delta
            wend = window.end(nsamps)*self.delta
            ax0.axvline(wbeg,linewidth=2,color='r')
            ax0.axvline(wend,linewidth=2,color='r')
            # windowed data
            d2 = self.chop(window,copy=True)
            ax1 = plt.subplot(gs[1])
            ax1.plot(d2.t()+wbeg,d2.data[0])
            ax1.plot(d2.t()+wbeg,d2.data[1])
            # particle  motion
            lim = abs(d2.data.max()) * 1.1
            ax2 = plt.subplot(gs[2])
            ax2.axis('equal')
            ax2.plot(d2.data[1],d2.data[0])
            ax2.set_xlim([-lim,lim])
            ax2.set_ylim([-lim,lim])
            ax2.axes.get_xaxis().set_visible(False)
            ax2.axes.get_yaxis().set_visible(False)

        # show
        plt.show()
    
    def split(self,degrees,tlag,copy=False):
        """
        Applies splitting operator (phi,dt) to Pair.
        
        Rotates data so that trace1 is lined up with degrees (and trace2 90 degrees clockwise).
        Applies a relative time shift by the nearest even number of samples to tlag,
        trace1 is shifted tlag/2 forward in time, and trace2 tlag/2 backward in time.
        Then undoes the original rotation.
        """
        # convert time shift to nsamples -- must be even
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        # find appropriate rotation angle
        rangle = degrees - self.angle
        # apply splitting
        if copy == False:
            self.data = core.split(self.data,rangle,nsamps)
        else:
            dupe = self.copy()
            dupe.data = core.split(self.data,rangle,nsamps)
            return dupe
    
    def unsplit(self,degrees,tlag,copy=False):
        """
        Applies reverse splitting operator (phi,dt) to Pair.
        
        Rotates data so that trace1 is lined up with degrees (and trace2 90 degrees clockwise).
        Applies a relative time shift by the nearest even number of samples to tlag,
        trace1 is shifted tlag/2 backward in time, and trace2 tlag/2 forward in time.
        Then undoes the original rotation.
        """
        # convert time shift to nsamples -- must be even
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        # find appropriate rotation angle
        rangle = degrees - self.angle
        if copy == False:
            self.data = core.unsplit(self.data,rangle,nsamps)
        else:
            dupe = self.copy()
            dupe.data = core.unsplit(self.data,rangle,nsamps)
            return dupe
            
    def rotateto(self,degrees,copy=False):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        # find appropriate rotation angle
        rangle = degrees - self.angle        
        if copy == False:
            self.data = core.rotate(self.data,rangle)
            self.angle = degrees
        else:
            dupe = self.copy()
            dupe.data = core.rotate(self.data,rangle)
            dupe.angle = degrees
            return dupe
        

        
    def lag(self,tlag,copy=False):
        """
        Relative shift trace1 and trace2 by tlag seconds
        """
        # convert time shift to nsamples -- must be even
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        if copy == False:
            self.data = core.lag(self.data,nsamps)
        else:
            dupe = self.copy()
            dupe.data = core.lag(self.data,nsamps)
            return dupe
     
    def chop(self,window,copy=False):
        """
        Chop data around window
        """
        nsamps = window.width
        centre = self.centre() + window.offset 
        tukey = window.tukey
        # action
        if copy == False:
            self.data = core.chop(self.data,nsamps,centre,tukey)
        else:
            dupe = self.copy()
            dupe.data = core.chop(self.data,nsamps,centre,tukey)
            return dupe
        
    def window(self,time_centre,time_width,tukey=None):
        """
        Return a window object about time_centre with time_width.
        """
        tcs = int(time_centre / self.delta)
        offset = tcs - self.centre()
        # convert time to nsamples -- must be odd
        width = int(time_width / self.delta)
        width = width if width%2==1 else width + 1        
        return window.Window(width,offset,tukey=tukey)
        
    # def autowindow(self,time_centre=None):
    #     """
    #     Makes a guess based on energy near *time_centre* about a suitable window
    #
    #     *time centre* should be the shear wave pick at the centre of the energy packet.
    #     By default will use centre sample.
    #     """
    #     if time_centre is None:
    #         t0 = self.centre()
    #     else:
    #         t0 = int(time_centre / self.delta)
            
        
               
    def copy(self):
        return copy.copy(self)
        
    # def grideigval(self, maxshift=None, window=None, stepang=None, stepshift=None):
    #     """
    #     Return an EigenM (after Silver and Chan, 1991).
    #
    #     Uses the modified method for calculating degrees of freedom of Walsh et al. 2014.
    #     """
    #
    #     return EigenM(self)

