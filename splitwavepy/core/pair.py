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
        - cmplabels = ['cmp1','cmp2'] | list of strings
        - units = 's' (for labelling) | string
    
    Geometry Keyword Arguments (if in doubt don't use):
        - geom = 'geo' (N,E) / 'ray' (SV,SH) / 'cart' (X,Y)
        - cmpvecs = np.eye(2)
        - rcvloc = (lat,lon,r) / (x,y,z)
        - srcloc =  (lat,lon,r) / (x,y,z)
        - rayloc =  rcvloc # arbirary point specified by user
        - rayvec = [0,0,1] # wave plane normal

    Methods:
        # display
            - p()
            - p_pm() # plot particle motion
            - p_tr() # plot waveform trace
        # splitting
            - split(fast,lag)
            - unsplit(fast,lag)
        # useful
            - t()           # time
            - data()        # both x and y in one array
            - chop()        # windowed data
            - rotateto()
            - pca()         # principal component analysis
            - suggest_labels()
        # windowing
            - set_window(start,end,Tukey=None)
            
        # io
            - copy()
            - save()
        # in development
            - plot(pickMode) # pickmode in development
    """
    def __init__(self,*args,**kwargs):
        
        # important to do first
        self.delta = 1.
        if ('delta' in kwargs): self.delta = kwargs['delta']
               
                         
        # if no args make synthetic
        if len(args) == 0: 
            self.x, self.y = _synth(**kwargs)               
        # otherwise read in data                
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
                   
        # Pair must have a window
        self.set_window(**kwargs)
         
        # add geometry info 
        self.geom = 'geo'
        if ('geom' in kwargs): self.geom = kwargs['geom']
           
        self.cmpvecs = np.eye(2)  
        if ('cmpvecs' in kwargs): self.cmpvectors = kwargs['cmpvecs']
        
        self.rayvec = [0,0,1] # normal to shear plane, along Z-axis
        if ('rayvec' in kwargs): self.rayvec = kwargs['rayvec']
        
        # source and receiver location info
        if ('srcloc' in kwargs): self.srcloc = kwargs['srcloc']     
        if ('rcvloc' in kwargs): self.rcvloc = kwargs['rcvloc']
        if ('rayloc' in kwargs): self.raylic = kwargs['rayloc']

        # labels
        self.units = 's'   
        if ('units' in kwargs): self.units = kwargs['units']      
        self.cmplabels = self.set_labels()
        if ('cmplabels' in kwargs): self.cmplabels = kwargs['cmplabels']   
        # A user defined name # maybe useful?
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
    
    def set_labels(self):
        if np.allclose(self.cmpvecs,np.eye(2)):
            if self.geom == 'geo': self.cmplabels = ['North','East']
            if self.geom == 'ray': self.cmplabels = ['Vertical','Horizontal']
            if self.geom == 'cart': self.cmplabels = ['X','Y']
            return
        # if reached here use the default label
        lab1,lab2 = self.cmpangs()
        self.cmplabels = ['Cmp '+str(round(lab1)),'Cmp ' +str(round(lab2))]
        
        
    def split(self,fast,lag):
        """
        Applies splitting operator.
        
        .. warning:: shortens trace length by *lag*.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag,self.delta,mode='even')
        # find appropriate rotation angle
        origangs=self.cmpangs()
        self.rotateto(0)
        # apply splitting
        self.x, self.y = core.split(self.x,self.y,fast,samps)
        self.rotateto(origangs[0])
        
        # ang = np.deg2rad(fast)
        # fastcmps = np.array([[np.cos(ang),-np.sin(ang)],
        #                      [np.sin(ang), np.cos(ang)]])
        # rot = np.dot(fastcmps,self.cmpvecs.T)
        # rangle = np.rad2deg(np.arccos(rot[0,0]))
        # print(rangle,samps)
        # apply splitting
        # self.x, self.y = core.split(self.x,self.y,rangle,samps)
           
    def unsplit(self,fast,lag):
        """
        Reverses splitting operator.
        
        .. warning:: shortens trace length by *lag*.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag,self.delta,mode='even')
        # find appropriate rotation angle
        origangs=self.cmpangs()
        self.rotateto(0)
        # apply splitting
        self.x, self.y = core.unsplit(self.x,self.y,fast,samps)
        self.rotateto(origangs[0])
        
        # ang = np.deg2rad(fast)
        # fastcmps = np.array([[np.cos(ang),-np.sin(ang)],
        #                      [np.sin(ang), np.cos(ang)]])
        # rot = np.dot(fastcmps,self.cmpvecs.T)
        # rangle = np.rad2deg(np.arccos(rot[0,0]))
        # print(rangle,samps)
        # self.x, self.y = core.unsplit(self.x,self.y,rangle,samps)
        
    def rotateto(self,degrees):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        # find appropriate rotation matrix
        ang = np.deg2rad(degrees)
        # define the new cmpvecs
        rotcmpvecs = np.array([[np.cos(ang),-np.sin(ang)],
                            [np.sin(ang), np.cos(ang)]])
        # find the rotation matrix. 
        # Linear algebra: if a and b are rotation matrices, 
        # then: a.T = inv(a) and b.T = inv(b)
        # then: dot(a.T,a) = dot(b.T,b) = I
        # and (multiply by b): dot(dot(b,a.T),a) = b.
        # i.e., dot(b,a.T) is the rotation matrix that converts a to b.
        rot = np.dot(rotcmpvecs,self.cmpvecs.T)
        # rotate data to suit
        xy = np.dot(rot,self.data())
        self.x, self.y = xy[0],xy[1]
        self.cmpvecs = rotcmpvecs
        self.set_labels()
        
        
    def pol(self):
        """Return principal component orientation"""
        # rotate to zero
        rot = self.cmpvecs.T
        xy = np.dot(rot,self.chop())
        _,eigvecs = core.eigcov(xy)
        x,y = eigvecs[:,0]
        pol = np.rad2deg(np.arctan2(y,x))
        return pol

        
    def eigen(self,window=None):
        self.eigvals, self.eigvecs = core.eigcov(self.data())


        

    def cmpangs(self):
        cmp1 = self.cmpvecs[:,0]
        cmp2 = self.cmpvecs[:,1]
        def getang(c) : return np.rad2deg(np.arctan2(c[1],c[0]))
        return getang(cmp1),getang(cmp2)
        
    

# Geometry stuff


# Windowing
                
    def set_window(self,*args,**kwargs):
        """
        Return a window object at user defined start and end times.
        
        args
        - Window
        - start,end
        
        kwargs
        - tukey
        
        The window will be adjusted to ensure it occupies an odd number 
        of samples.
        """
                
        # if Window provided
        if 'window' in kwargs:  
            if isinstance(kwargs['window'],Window):
                self.window = kwargs['window']
                return
            else:
                raise TypeError('expecting a window')
        
        # if no arguments provided
        if len(args) == 0:
            width = self.nsamps() / 3
            self.window = Window(width)
            return
        # if start/end given
        if len(args) == 2:
            start, end = args 
            time_centre = (start + end)/2
            time_width = end - start
            tcs = core.time2samps(time_centre)
            offset = tcs - self.centre()
            # convert time to nsamples -- must be odd
            width = core.time2samps(time_width,delta,'odd')       
            self.window = Window(width,offset,**kwargs) 
            return
        else:
            raise Exception ('unexpected number of arguments')     
    
    def chop(self):
        """
        Chop data to window
        """
        x,y = core.chop(self.x,self.y,window=self.window)
        return x,y
        
    def chopt(self):
        """
        Chop time to window
        """        
        t = core.chop(self.t(),window=self.window)
        return t

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


    # useful?
    def nsamps(self):
        return self.x.size
    
    def centresamp(self):
        return int(self.x.size/2)
        
    def centretime(self):
        return int(self.x.size/2) * self.delta
        
    def power(self):
        return self.x**2+self.y**2
                


        

            
        
               
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
    
    def p_tr(self,**kwargs):
        """Plot traces"""
        ax = plot.trace(self.x,self.y,time=self.t(),**kwargs)
        plt.show()
    
    def p_pm(self,**kwargs):
        """Plot particle motion"""
        ax = plot.particle(self.x,self.y,cmplabels=self.cmplabels,**kwargs)
        plt.show()
           
    def plot(self,**kwargs):
        """
        Plot trace data and particle motion
        """
        from matplotlib import gridspec
        fig = plt.figure(figsize=(12, 3)) 
               
        if 'window' in kwargs and kwargs['window'] is True:                        
            gs = gridspec.GridSpec(1, 3, width_ratios=[3,1,1])
            # trace with window markers
            ax0 = plt.subplot(gs[0])
            plot.trace(self.x,self.y,time=self.t(),window=self.window,ax=ax0)
            # windowed data
            ax1 = plt.subplot(gs[1])
            chopx, chopy = self.chop()
            chopt = self.chopt()
            plot.trace( chopx, chopy, time=chopt, ax=ax1)
            # particle  motion
            ax2 = plt.subplot(gs[2])
            # ax2.axis('equal')
            plot.particle( chopx, chopy, ax=ax2)                        
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
    noise = 0.001
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
    


