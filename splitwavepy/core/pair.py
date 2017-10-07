# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core
from ..plotting import plot
from .window import Window
from ..core import io

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
        - cmplabels = ['Comp1','Comp2'] | list of strings
        - units = 's' (for labelling) | string
    
    Geometry Keyword Arguments (if in doubt don't use):
        - geom = 'geo' (N,E) / 'ray' (SV,SH) / 'cart' (X,Y)
        - cmpvecs = np.eye(2)
        - rcvloc = (lat,lon,r) / (x,y,z)
        - srcloc =  (lat,lon,r) / (x,y,z)
        - rayloc =  rcvloc # must be somewhere on a sensible raypath
        - rayvector = [0,0,1] 

    Methods:
        - plot()
        - plot(pickMode)
    """
    def __init__(self,*args,**kwargs):
        
        # important to do first
        self.delta = 1.
        if ('delta' in kwargs): self.delta = kwargs['delta']        
        if ('window' in kwargs): self.window = kwargs['window']
                        
        # make synthetic
        if len(args) == 0:
            # generate                                                       
            self.x, self.y = _synth(**kwargs)  
             
        # read in data                
        elif len(args) == 2:
            if not (isinstance(args[0],np.ndarray) & isinstance(args[0],np.ndarray)):
                raise TypeError('expecting numpy arrays')         
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
           
        self.cmpvecs = np.eye(2)  
        if ('cmpvecs' in kwargs): self.cmpvectors = kwargs['cmpvecs']
        
        self.rayvector = [0,0,1] # normal to shear plane, along Z-axis
        
        if ('srcloc' in kwargs): self.srcloc = kwargs['srcloc']     
        if ('rcvloc' in kwargs): self.rcvloc = kwargs['rcvloc']
        if ('rayloc' in kwargs): self.raylic = kwargs['rayloc']

        # labels
        self.units = 's'   
        if ('units' in kwargs): self.units = kwargs['units']      
        self.cmplabels = ['Comp1','Comp2']
        if ('cmplabels' in kwargs): self.cmplabels = kwargs['cmplabels']   
        # A user defined name
        self.name = 'untitled'
        if ('name' in kwargs): self.name = kwargs['name']    
        # Backup the command used to produce this object
        self.args = args
        self.kwargs = kwargs

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
    
    def split(self,degrees,lag):
        """
        Applies splitting operator (phi,dt) to Pair.
        
        Rotates data so that trace1 is lined up with degrees (and trace2 90 degrees clockwise).
        Applies a relative time shift by the nearest even number of samples to tlag,
        trace1 is shifted tlag/2 forward in time, and trace2 tlag/2 backward in time.
        Then undoes the original rotation.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag,self.delta,mode='even')
        # find appropriate rotation angle
        rangle = degrees - self.cmpangles()[0]
        # apply splitting
        self.x, self.y = core.split(self.x,self.y,rangle,nsamps)
           
    def unsplit(self,degrees,lag):
        """
        Applies reverse splitting operator (phi,dt) to Pair.
        
        Rotates data so that trace1 is lined up with degrees (and trace2 90 degrees clockwise).
        Applies a relative time shift by the nearest even number of samples to tlag,
        trace1 is shifted tlag/2 backward in time, and trace2 tlag/2 forward in time.
        Then undoes the original rotation.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag,self.delta,mode='even')
        # find appropriate rotation angle
        rangle = degrees - self.cmpangs()[0]
        self.x, self.y = core.unsplit(self.x,self.y,rangle,nsamps)
        
    def rotateto(self,degrees):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        # find appropriate rotation matrix
        ang = np.deg2rad(degrees)
        # define the new cmpvecs
        newcmpvecs = np.array([[np.cos(ang),-np.sin(ang)],
                        [np.sin(ang), np.cos(ang)]])
        # find the rotation matrix. 
        # Linear algebra: if a and b are rotation matrices, 
        # which have the useful property: a.T = inv(a)
        # then: dot(a.T,a) = dot(b.T,b) = I
        # and (multiply by b): dot(dot(b,a.T),a) = b.
        # i.e., dot(b,a.T) is the rotation matrix that converts a to b.
        rot = np.dot(newcmpvecs,self.cmpvecs.T)
        # rotate data to suit
        xy = np.dot(rot,self.data())
        self.x, self.y = xy[0],xy[1]
        self.cmpvecs = newcmpvecs

                
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
        
    def cmpangs(self):
        cmp1 = self.cmpvecs[:,0]
        cmp2 = self.cmpvecs[:,1]
        def getang(c) : return np.rad2deg(np.arctan2(c[1],c[0]))
        return getang(cmp1),getang(cmp2)

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
           
    def plot(self,**kwargs):
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
        #
        # # plot mode
        # if args[0] == 'pickmode':
        #     windowpicker = plot.WindowPicker(fig,ax0,**kwargs)
        #     windowpicker.connect()
        #     # self.getWindow()
                     
        # show
        plt.tight_layout()
        plt.show()      
        
def _synth(**kwargs):
    """return ricker wavelet synthetic data"""
    
    # defaults
    pol = 0.
    delta = 1.
    fast = 0.
    lag = 0.
    noise = 0.03
    nsamps = 501
    width = 16.
    
    # override defaults
    if ('pol' in kwargs): pol = kwargs['pol']   
    if ('delta' in kwargs): delta = kwargs['delta']  
    if ('fast' in kwargs): fast = kwargs['fast']
    if ('lag' in kwargs): lag = kwargs['lag']   
    if ('noise' in kwargs): noise = kwargs['noise']   
    if ('nsamps' in kwargs): nsamps = kwargs['nsamps']   
    if ('width' in kwargs): width = kwargs['width']   

    nsamps = int(nsamps)  
    x = signal.ricker(nsamps, width) + core.noise(nsamps,noise,int(width/4))
    y = core.noise(nsamps,noise,width/4)    
    # rotate to polarisation 
    # negative because we are doing the active rotation of data, whereas
    #core.rotate does the passive transormation of the co-ordinate system
    # (generally more common in splitting).
    x,y = core.rotate(x,y,-pol)    
    # add any splitting -- lag samples must be even
    slag = core.time2samps(lag,delta,mode='even')
    x,y = core.split(x,y,fast,slag)
    
    return x,y
    


