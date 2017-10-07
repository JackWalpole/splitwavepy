# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core
from ..plotting import plot
from .window import Window

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import copy

class Pair:
    """
    The Pair: work with 2-component data.
        
    Usage: Pair()     => create Pair of synthetic data
           Pair(x,y) => creates Pair from two traces stored in numpy arrays x and y.
    
    Keyword Arguments:
        - delta = 1. (sample interval) [default] | float
    
    Naming Keyword Arguments:
        - name = 'untitled' (should be unique identifier) | string
        - units = 's' (for labelling) | string
    
    Geometry Keyword Arguments (if in doubt don't use):
        - geom = 'geo' (N,E) / 'ray' (SV,SH) / 'cart' (X,Y)
        - cmpvectors = np.eye(2)
        - rcvloc = (lat,lon,r) / (x,y,z)
        - srcloc =  (lat,lon,r) / (x,y,z)
        - rayloc =  rcvloc # must be somewhere on a sensible raypath
        - rayvector = [0,0,1] 

    Methods:
        - plot()
        - plot(pickMode)
    """
    def __init__(self,*args,**kwargs):
        
        # defaults
        self.delta = 1.
        if ('delta' in kwargs): self.delta = kwargs['delta']
        
        # labels
        self.units = 's'   
        if ('units' in kwargs): self.units = kwargs['units']
        
        # A user defined name
        kwargs['name'] = 'untitled'
        if ('name' in kwargs): self.name = kwargs['name']
        
        # make synthetic
        if len(args) == 0:
            # convert time2samps -- must be even
            if ('lag' in kwargs): 
                kwargs['lag'] = core.time2samps(kwargs['lag'])
            # generate                                                       
            self.x, self.y = _synth(**kwargs)   
        # read in data                
        elif len(args) == 2:            
            self.x, self.y = args[0], args[1]
        else: 
            raise Exception('Unexpected number of arguments')
                    
        # some sanity checks
        if self.x.ndim != 1:
            raise Exception('data must be one dimensional')
        if self.x.size%2 == 0:
            raise Exception('data must have odd number of samples')
        if (self.x.size != self.y.size):
            raise Exception('x and y must be the same length')
        
        # add geometry info 
        self.geom = 'geo'
        if ('geom' in kwargs): self.geom = kwargs['geom']
           
        self.cmpvectors = np.eye(2)  
        if ('cmpvectors' in kwargs): self.cmpvectors = kwargs['cmpvectors']
        
        self.rayvector = [0,0,1] # normal to shear plane, along Z-axis
        
        if ('srcloc' in kwargs): self.srcloc = kwargs['srcloc']     
        if ('rcvloc' in kwargs): self.rcvloc = kwargs['rcvloc']
        self.rayloc = self.rcvloc
        if ('rayloc' in kwargs): self.raylic = kwargs['rayloc']


    # methods

    def t(self):
        return np.arange(self.x.size) * self.delta
  
    def data(self):
        return np.vstack((self.x,self.y))
    
    def get_labels(self):
        if self.geom == 'geo': return ['North','East']
        if self.geom == 'ray': return ['Vertical','Horizontal']
        if self.geom == 'cart': return ['X','Y']
    
    def pca(self):
        """Return orientation of principal component"""
        return core.pca(self.x,self.y)
    
    def split(self,degrees,tlag):
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
        self.x, self.y = core.split(self.x,self.y,rangle,nsamps)
           
    def unsplit(self,degrees,tlag):
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
        self.x, self.y = core.unsplit(self.x,self.y,rangle,nsamps)
        
    def rotateto(self,degrees):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        # find appropriate rotation angle
        rangle = degrees - self.angle 
        self.x, self.y = core.rotate(self.x,self.y,rangle)
        self.angle = degrees
                    

    def setWindow(self,start,end,tukey=None):
        """
        Return a window object at user defined start and end times.
        
        The window will be adjusted to ensure it occupies an odd number 
        of samples.
        """
        time_centre = (start + end)/2
        time_width = end - start
        tcs = int(time_centre / self.delta)
        offset = tcs - self.centre()
        # convert time to nsamples -- must be odd
        width = int(time_width / self.delta)
        width = width if width%2==1 else width + 1        
        self.window = Window(width,offset,tukey=tukey)        
     
    def chop(self,window):
        """
        Chop data around window
        """
        self.x, self.y = core.chop(self.x,self.y,window=window)
        


    # useful?
    def nsamps(self):
        return self.x.size
    
    def centresamp(self):
        return int(self.x.size/2)
        
    def centretime(self):
        return int(self.x.size/2) * self.delta
        
    def power(self):
        return self.x**2+self.y**2
                


        
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
        return copy.deepcopy(self)
        
        
    # Comparison
    
    def __eq__(self, other) :
        # check same class
        if self.__class__ != other.__class__: return False
        # check same keys
        if self.__dict__.keys() != other.__dict__.keys(): return False
        # check same values
        for key in self.__dict__.keys():
            if np.all( self.__dict__[key] != other.__dict__[key]): return False
        # if reached here then the same
        return True
    
    # plotting
    
    def pt(self,**kwargs):
        """Plot traces"""
        ax = plot.trace(self.x,self.y,time=self.t(),**kwargs)
        plt.show()
    
    def ppm(self,**kwargs):
        """Plot particle motion"""
        ax = plot.particle(self.x,self.y,**kwargs)
        plt.show()
           
    def plot(self,*args,**kwargs):
        """
        Plot trace data and particle motion
        """
        from matplotlib import gridspec
        fig = plt.figure(figsize=(12, 3)) 
               
        if 'window' in kwargs:                        
            if kwargs['window'] is True:
                kwargs['window'] = self.window
            window = kwargs['window']
            gs = gridspec.GridSpec(1, 3, width_ratios=[3,1,1])
            # trace with window markers
            ax0 = plt.subplot(gs[0])
            plot.trace(self.x,self.y,time=self.t(),window=window,ax=ax0)
            # windowed data
            nsamps = self.nsamps()
            wbeg = window.start(nsamps)*self.delta
            d2 = self.copy()
            d2.chop(window)
            ax1 = plt.subplot(gs[1])
            plot.trace(d2.x,d2.y,time=d2.t()+wbeg,ax=ax1)
            # particle  motion
            ax2 = plt.subplot(gs[2])
            # ax2.axis('equal')
            plot.particle(d2.x,d2.y,ax=ax2)                        
        else:
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
            # trace
            ax0 = plt.subplot(gs[0])
            plot.trace(self.x,self.y,time=self.t(),ax=ax0)
            # particle  motion
            ax1 = plt.subplot(gs[1])
            plot.particle(self.x,self.y,ax=ax1)
            
        # plot mode
        if args[0] == 'pickmode':            
            windowpicker = plot.WindowPicker(fig,ax0,**kwargs)
            windowpicker.connect()
            # self.getWindow()
                     
        # show
        plt.tight_layout()
        plt.show()      
        
def _synth(**kwargs):
    """return ricker wavelet synthetic data"""
    
    if ('pol' in kwargs):
        pol = kwargs['pol']
    else:
        pol = 0.
        
    if ('fast' in kwargs):
        fast = kwargs['fast']
    else:
        fast = 0.
        
    if ('lag' in kwargs):
        lag = kwargs['lag']
    else:
        lag = 0
        
    if ('noise' in kwargs):
        noise = kwargs['noise']
    else:
        noise = 0.03
        
    if ('nsamps' in kwargs):
        nsamps = kwargs['nsamps']
    else:
        nsamps = 501
        
    if ('width' in kwargs):
        width = kwargs['width']
    else:
        width = 16.
        
    if ('window' in kwargs):
        window = kwargs['window']
    else:
        window = Window(width*3)

    nsamps = int(nsamps)
    
    x = signal.ricker(nsamps, width) + core.noise(nsamps,noise,width/4)
    y = core.noise(nsamps,noise,width/4)
    
    # rotate to polarisation
    x,y = core.rotate(x,y,-pol)
    
    # add any splitting -- this will reduce nsamps
    x,y = core.split(x,y,fast,lag)
    
    return x,y
    


