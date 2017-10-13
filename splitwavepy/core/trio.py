# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, core3d, geom, io
from .pair import Pair
from .window import Window

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec

class Trio:
    """
    The Trio: work with 3-component data.
        
    Usage: Trio()     => create Trio of synthetic data
           Trio(data) => creates Trio from two traces stored as rows in numpy array data
           Trio(x,y,z) => creates Trio from two traces stored in numpy arrays x and y.
    
    Keyword Arguments:
        - delta = 1. (sample interval) [default] | float
        - units = 's' (for labelling) | string
        - geom = 'cart' (x,y,z) [default] | 'geo' (az,inc,r) | 'ray' (P,SH,SV) 
        - window = None (default) | Window object
        - angle
    
    Advanced Keyword Arguments (if in doubt don't use):
        - xyz = np.ones(3) | custom numpy array
        - rcvloc = None
        - srcloc = None
    """
    def __init__(self,*args,**kwargs):
        
        # important to do first
        self.delta = 1.
        if ('delta' in kwargs): self.delta = kwargs['delta']
                       
        # if no args make synthetic
        if len(args) == 0: 
            self.x, self.y, self.z = _synth(**kwargs)               
        # otherwise read in data                
        elif len(args) == 3:
            if not (isinstance(args[0], np.ndarray) & 
                    isinstance(args[1], np.ndarray) &
                    isinstance(args[2], np.ndarray)):
                raise TypeError('expecting numpy arrays')         
            self.x, self.y, self.z = args[0], args[1], args[2]
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
           
        self.cmpvecs = np.eye(3)  
        if ('cmpvecs' in kwargs): self.cmpvecs = kwargs['cmpvecs']
        
        self.rayvec = [0,0,1] # normal to shear plane, along Z-axis
        if ('rayvec' in kwargs): self.rayvec = kwargs['rayvec']
        
        # source and receiver location info
        if ('srcloc' in kwargs): self.srcloc = kwargs['srcloc']     
        if ('rcvloc' in kwargs): self.rcvloc = kwargs['rcvloc']
        if ('rayloc' in kwargs): self.raylic = kwargs['rayloc']

        # labels
        self.units = 's'
        self.set_labels()
        self.name = 'untitled'
        if ('units' in kwargs): self.units = kwargs['units']      
        if ('cmplabels' in kwargs): self.cmplabels = kwargs['cmplabels']
        if ('name' in kwargs): self.name = kwargs['name']    
        # Backup the command used to produce this object
        self.args = args
        self.kwargs = kwargs

    # METHODS
        
    # def split(self,fast,lag):
    #     """
    #     Applies splitting operator.
    #
    #     .. warning:: shortens trace length by *lag*.
    #     """
    #     # convert time shift to nsamples -- must be even
    #     samps = core.time2samps(lag,self.delta,mode='even')
    #     # find appropriate rotation angle
    #     origangs=self.cmpangs()
    #     self.rotateto(0)
    #     # apply splitting
    #     self.x, self.y = core.split(self.x,self.y,fast,samps)
    #     self.rotateto(origangs[0])
    #
    # def unsplit(self,fast,lag):
    #     """
    #     Reverses splitting operator.
    #
    #     .. warning:: shortens trace length by *lag*.
    #     """
    #     # convert time shift to nsamples -- must be even
    #     samps = core.time2samps(lag,self.delta,mode='even')
    #     # find appropriate rotation angle
    #     origangs=self.cmpangs()
    #     self.rotateto(0)
    #     # apply splitting
    #     self.x, self.y = core.unsplit(self.x,self.y,fast,samps)
    #     self.rotateto(origangs[0])
        
    def rotateto(self,vecs):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        # define the new cmpvecs
        backoff = self.cmpvecs.T
        self.cmpvecs = rotcmpvecs
        # find the rotation matrix. 
        # Linear algebra: if a and b are rotation matrices, 
        # then: a.T = inv(a) and b.T = inv(b)
        # then: dot(a.T,a) = dot(b.T,b) = I
        # and (multiply by b): dot(dot(b,a.T),a) = b.
        # i.e., dot(b,a.T) is the rotation matrix that converts a to b.
        rot = np.dot(self.cmpvecs,backoff)
        # rotate data to suit
        xyz = np.dot(rot,self.data())
        self.x, self.y,self.z = xyz[0],xyz[1],xyz[2]
        # reset label
        # if reached here use the default label
        # lab1,lab2,lab3 = self.cmpangs()
        # self.set_labels()

    def rotz(self,degs):
        """Rotate about z axis."""
        rads = math.radians(degs)
        rotateto( geom.rz( self.cmpvecs,rads))
        
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
            width = core.odd(self._nsamps() / 3)
            self.window = Window(width)
            return
        # if start/end given
        if len(args) == 2:
            start, end = args 
            time_centre = (start + end)/2
            time_width = end - start
            tcs = core.time2samps(time_centre,self.delta)
            offset = tcs - self._centresamp()
            # convert time to nsamples -- must be odd
            width = core.time2samps(time_width,self.delta,'odd')       
            self.window = Window(width,offset,**kwargs) 
            return
        else:
            raise Exception ('unexpected number of arguments')

    def set_labels(self,*args):
        if len(args) == 0:
            if np.allclose(self.cmpvecs,np.eye(3),atol=1e-02):
                if self.geom == 'geo': self.cmplabels = ['North','East','Z'] #
                elif self.geom == 'ray': self.cmplabels = ['Vertical','Horizontal','Ray']
                elif self.geom == 'cart': self.cmplabels = ['X','Y','Z']
                else: self.cmplabels = ['Comp1','Comp2','Comp3']
                return
            # if reached here we have a non-standard orientation
            # a1,a2 = self.cmpangs()
            # lbl1 = str(round(a1))+r' ($^\circ$)'
            # lbl2 = str(round(a2))+r' ($^\circ$)'
            # self.cmplabels = [lbl1,lbl2]
            self.cmplabels = ['Comp1','Comp2','Comp3']
            return
        elif len(args) == 1:
            if not isinstance(args[0],list): raise TypeError('expecting a list')
            if not len(args[0]) == 3: raise Exception('list must be length 3')
            if not (isinstance(args[0][0],str) and 
                    isinstance(args[0][1],str) and
                    isinstance(args[0][2],str)):
                raise TypeError('cmplabels must be a list of strings')
            self.cmplabels = args[0]
            return
        else:
            raise Exception('unexpected number of arguments')



                

        
    # Utility 
    
    def t(self):
        return np.arange(self.x.size) * self.delta
  
    def data(self):
        return np.vstack(( self.x, self.y, self.z))
        
    def pol(self):
        """Return principal component orientation"""
        # rotate to zero
        rot = self.cmpvecs.T
        data = self.chop().data()
        xy = np.dot(rot,data)
        _,eigvecs = core.eigcov(xy)
        pol = eigvecs[:,0]
        return pol

    def eigdata(self):
        """Return to maximum, intermediate, and minimum directions."""
        # rotate to I
        rot = self.cmpvecs.T
        data = self.chop().data()
        xyz = np.dot(rot,data)
        _,eigvecs = core.eigcov(xyz)
        rotateto(self,eigvecs)
        
    def eigvals(self):
        """Return principal component vector."""
        # rotate to I
        rot = self.cmpvecs.T
        data = self.chop().data()
        xyz = np.dot(rot,data)
        eigvals,_ = core.eigcov(xyz)
        return(eigvals)  
        

        

    def power(self):
        return self.x**2,self.y**2,self.z**2

    # def cmpangs(self):
    #     cmp1 = self.cmpvecs[:,0]
    #     cmp2 = self.cmpvecs[:,1]
    #     def getang(c) : return np.rad2deg(np.arctan2(c[1],c[0]))
    #     return getang(cmp1),getang(cmp2)        
    
    def chop(self):
        """
        Chop data to window
        """
        chop = self.copy()
        chop.x, chop.y, chop.z = core.chop(self.x,self.y,self.z,window=self.window)
        self.window.offset = 0
        return chop
        
    def chopt(self):
        """
        Chop time to window
        """        
        t = core.chop(self.t(),window=self.window)
        return t
    
    def wbeg(self):
        """
        Window start time.
        """
        sbeg = self.window.start(self._nsamps())
        return sbeg * self.delta
    
    def wend(self):
        """
        Window end time.
        """
        send = self.window.end(self._nsamps())
        return send * self.delta
        
    def wwidth(self):
        """
        Window width.
        """
        return self.window.width * self.delta

    # Plotting
 
 # def plot(self,window=None):
 #     """
 #     Plot trace data and particle motion
 #     """
 #     from matplotlib import gridspec
 #     fig = plt.figure(figsize=(12, 3))
 #     if window is None:
 #         gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
 #         # trace
 #         ax0 = plt.subplot(gs[0])
 #         plot.trace(self.x,self.y,self.z,time=self.t(),ax=ax0)
 #         # particle  motion
 #         ax1 = plt.subplot(gs[1],projection='3d')
 #         plot.particle(self.x,self.y,self.z,ax=ax1)
 #     else:
 #         gs = gridspec.GridSpec(1, 3, width_ratios=[3,1,1])
 #         # trace with window markers
 #         ax0 = plt.subplot(gs[0])
 #         plot.trace(self.x,self.y,self.z,time=self.t(),window=window,ax=ax0)
 #         # windowed data
 #         d2 = self.copy()
 #         d2.chop(window)
 #         ax1 = plt.subplot(gs[1])
 #         nsamps = self.nsamps()
 #         wbeg = window.start(nsamps)*self.delta
 #         plot.trace(d2.x,d2.y,d2.z,time=d2.t()+wbeg,ax=ax1)
 #         # particle  motion
 #         ax2 = plt.subplot(gs[2],projection='3d')
 #         plot.particle(d2.x,d2.y,d2.z,ax=ax2)
 #     # show
 #     plt.tight_layout()
 #     plt.show()
 #
 #
 # def pt(self,**kwargs):
 #     """Plot traces"""
 #     ax = plot.trace(self.x,self.y,self.z,time=self.t(),**kwargs)
 #     plt.show()
 #
 # def ppm(self,**kwargs):
 #     """Plot particle motion"""
 #     ax = plot.particle(self.x,self.y,self.z,**kwargs)
 #     plt.show()
     
                
    def plot(self,**kwargs):
        """
        Plot trace data and particle motion
        """

        fig = plt.figure(figsize=(12, 3))     
        gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
        
        # trace
        ax0 = plt.subplot(gs[0])
        self._ptr( ax0, **kwargs)
        
        # particle  motion
        ax1 = plt.subplot(gs[1],projection='3d')
        self._ppm( ax1, **kwargs)   
                                 
        # show
        plt.tight_layout()
        plt.show()      

    def _ptr( self, ax, **kwargs):
        """Plot trace data on *ax* matplotlib axis object.
        """    
        # plot data
        t = self.t()
        
        # set labels
        if 'cmplabels' not in kwargs: kwargs['cmplabels'] = self.cmplabels
        ax.plot( t, self.x, label=kwargs['cmplabels'][0])
        ax.plot( t, self.y, label=kwargs['cmplabels'][1])
        ax.plot( t, self.z, label=kwargs['cmplabels'][2])
        ax.legend()
    
        # set limits
        lim = np.abs(self.data()).max() * 1.1
        if 'ylim' not in kwargs: kwargs['ylim'] = [-lim,lim]
        ax.set_ylim(kwargs['ylim'])
        if 'xlim' in kwargs: ax.set_xlim(kwargs['xlim'])
    
        # set axis label
        if 'units' not in kwargs: kwargs['units'] = 's'            
        ax.set_xlabel('Time (' + kwargs['units'] +')')

        # plot window markers
        # if 'window' not in kwargs and kwargs['window'] != False:
        if self.window.width < self._nsamps():
            w1 = ax.axvline(self.wbeg(),linewidth=1,color='k')
            w2 = ax.axvline(self.wend(),linewidth=1,color='k')
            w1.aname = 'w1'
            w2.aname = 'w2'       
        return

    def _ppm(self,ax,**kwargs):
        """Plot particle motion on *ax* matplotlib axis object.
        """
                
        # plot data
        xyz = self.chop().data()
        x,y,z = xyz[0], xyz[1], xyz[2]
        ax.plot( x, y, z)
    
        # set limit
        lim = np.abs(self.data()).max() * 1.1
        if 'lims' not in kwargs: kwargs['lims'] = [-lim,lim] 
        ax.set_aspect('equal')
        ax.set_xlim(kwargs['lims'])
        ax.set_ylim(kwargs['lims'])
        ax.set_zlim(kwargs['lims'])
    
        # set labels
        if 'cmplabels' not in kwargs: kwargs['cmplabels'] = self.cmplabels
        ax.set_xlabel(kwargs['cmplabels'][0])
        ax.set_ylabel(kwargs['cmplabels'][1])
        ax.set_zlabel(kwargs['cmplabels'][2])
        
        # side panel data
        ax.plot(x, y,-lim, zdir='z', alpha=0.3, color='g')
        ax.plot(x, z, lim, zdir='y', alpha=0.3, color='g')
        ax.plot(y, z,-lim, zdir='x', alpha=0.3, color='g')
        
        # turn off tick annotation
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.axes.zaxis.set_ticklabels([])
        return

    # I/O stuff  

    def save(self,filename):
        """
        Save Measurement for future referral
        """
        io.save(self,filename)
                       
    def copy(self):
        return io.copy(self)

        
    # Geometry stuff

    # def geom_to_geo():
    # def geom_to_ray():
    # def geom_to   
        
    # Hidden
    
    def _nsamps(self):
        return self.x.size

    def _centresamp(self):
        return int(self.x.size/2)
    
    def _centretime(self):
        return int(self.x.size/2) * self.delta

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

# Exterior   

def _synth(**kwargs):
    """return ricker wavelet synthetic data"""
    
    # defaults
    pol = 0.
    delta = 1.
    # fast = 0.
    # lag = 0.
    split = []
    noise = 0.001
    nsamps = 1001
    width = 32.
    
    # override defaults
    if ('pol' in kwargs): pol = kwargs['pol']   
    if ('delta' in kwargs): delta = kwargs['delta']  
    # if ('fast' in kwargs): fast = kwargs['fast']
    # if ('lag' in kwargs): lag = kwargs['lag']
    if ('split') in kwargs: split = kwargs['split']
    if ('noise' in kwargs): noise = kwargs['noise']   
    if ('nsamps' in kwargs): nsamps = kwargs['nsamps']   
    if ('width' in kwargs): width = kwargs['width'] 
    noisewidth = width/4  
    if ('noisewidth' in kwargs): noisewidth = kwargs['noisewidth']

    # initiate data with clean ricker wavelet
    nsamps = int(nsamps)  
    x = signal.ricker(nsamps, width)
    y = np.zeros(nsamps)
    
    # rotate to polarisation 
    # negative because we are doing the active rotation of data, whereas
    # core.rotate does the passive transormation of the co-ordinate system
    x,y = core.rotate(x,y,-pol)

    if isinstance(split,tuple):
        fast, lag = split
        # add any splitting -- lag samples must be even
        slag = core.time2samps(lag,delta,mode='even')
        x,y = core.split(x,y,fast,slag)
    elif isinstance(split,list):        
        for parms in split:
            fast, lag = parms
            # add any splitting -- lag samples must be even
            slag = core.time2samps(lag,delta,mode='even')
            x,y = core.split(x,y,fast,slag)
    
    # add noise - do this last to avoid splitting the noise
    x = x + core.noise(x.size,noise,int(noisewidth))    
    y = y + core.noise(x.size,noise,int(noisewidth))
    z = core.noise(x.size,noise,int(noisewidth))

    return x,y,z
