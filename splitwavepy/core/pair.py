# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, geom, io
from .window import Window

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
    def __init__(self,*args,**kwargs):
        
        # important to do first
        self.delta = 1.
        if ('delta' in kwargs): self.delta = kwargs['delta']
                       
                       
        # if no args make synthetic
        if len(args) == 0: 
            self.x, self.y = _synth(**kwargs)               
        # otherwise read in data                
        elif len(args) == 2:
            if not (isinstance(args[0],np.ndarray) & isinstance(args[1],np.ndarray)):
                raise TypeError('expecting numpy arrays')         
            self.x, self.y = args[0], args[1]
        else: 
            raise Exception('Unexpected number of arguments')
                    
        # some sanity checks
        if self.x.ndim != 1: raise Exception('data must be one dimensional')
        if self.x.size%2 == 0: raise Exception('data must have odd number of samples')
        if (self.x.size != self.y.size): raise Exception('x and y must be the same length')
                   
        # Pair must have a window
        self.set_window(**kwargs)
         
        # add geometry info 
        self.geom = 'geo'
        if ('geom' in kwargs): self.geom = kwargs['geom']
           
        self.cmpvecs = np.eye(2)  
        if ('cmpvecs' in kwargs): self.cmpvecs = kwargs['cmpvecs']
        
        # if pol specified set
        if ('pol' in kwargs): 
            self.set_pol(kwargs['pol'])
        else:
            self.set_pol()
        
        # self.rayvec = [0,0,1] # normal to shear plane, along Z-axis
        # if ('rayvec' in kwargs): self.rayvec = kwargs['rayvec']
        # Always just assume ray vector is normal to components
        
        # source and receiver location info
        if ('srcloc' in kwargs): self.srcloc = kwargs['srcloc']     
        if ('rcvloc' in kwargs): self.rcvloc = kwargs['rcvloc']
        if ('rayloc' in kwargs): self.raylic = kwargs['rayloc']

        # labels
        self.units = 's'   
        if ('units' in kwargs): self.units = kwargs['units']      
        self.set_labels()
        if ('cmplabels' in kwargs): self.cmplabels = kwargs['cmplabels']
        # A user defined name # maybe useful?
        self.name = 'untitled'
        if ('name' in kwargs): self.name = kwargs['name']    
        # Backup the command used to produce this object
        self.args = args
        self.kwargs = kwargs

    # METHODS
      
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
       
    def rotateto(self,degrees):
        """
        Rotate traces so that cmp1 lines up with *degrees*
        """
        # find appropriate rotation matrix
        ang = math.radians(degrees)
        # define the new cmpvecs
        backoff = self.cmpvecs
        self.cmpvecs = np.array([[ math.cos(ang),-math.sin(ang)],
                                 [ math.sin(ang), math.cos(ang)]])
        rot = np.dot(self.cmpvecs.T, backoff)
        # rotate data
        xy = np.dot(rot, self.data())
        self.x, self.y = xy[0], xy[1]
        # reset label
        self.set_labels()
        
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
            self.window = self._construct_window(start,end,**kwargs) 
            return
        else:
            raise Exception ('unexpected number of arguments')

    def set_labels(self,*args):
        if len(args) == 0:
            if np.allclose(self.cmpvecs,np.eye(2),atol=1e-02):
                if self.geom == 'geo': self.cmplabels = ['North','East']
                elif self.geom == 'ray': self.cmplabels = ['Vertical','Horizontal']
                elif self.geom == 'cart': self.cmplabels = ['X','Y']
                else: self.cmplabels = ['Comp1','Comp2']
                return
            # if reached here we have a non-standard orientation
            a1,a2 = self.cmpangs()
            lbl1 = str(round(a1))+r' ($^\circ$)'
            lbl2 = str(round(a2))+r' ($^\circ$)'
            self.cmplabels = [lbl1,lbl2]
            return
        elif len(args) == 1:
            if not isinstance(args[0],list): raise TypeError('expecting a list')
            if not len(args[0]) == 2: raise Exception('list must be length 2')
            if not (isinstance(args[0][0],str) and isinstance(args[0][1],str)):
                raise TypeError('cmplabels must be a list of strings')
            self.cmplabels = args[0]
            return
        else:
            raise Exception('unexpected number of arguments')
            
    def set_pol(self,*args):
        if len(args) == 0:
            self.pol = self.get_pol()
        elif len(args) == 1:
            self.pol = float(args[0])
        else:
            raise Exception('Unexpected number of arguments')
        return
    
    # Utility 
    
    def t(self):
        return np.arange(self.x.size) * self.delta
  
    def data(self):
        return np.vstack((self.x,self.y))

    def get_pol(self):
        """Return principal component orientation"""
        # rotate to zero
        rot = self.cmpvecs.T
        data = self.chop().data()
        xy = np.dot(rot,data)
        _,eigvecs = core.eigcov(xy)
        x,y = eigvecs[:,0]
        pol = np.rad2deg(np.arctan2(y,x))
        return pol
        
    def eigen(self,window=None):
        self.eigvals, self.eigvecs = core.eigcov(self.data())
        
    def power(self):
        return self.x**2,self.y**2
        
    # def snrRH(self):
    #     data = self.copy()
    #     data.rotateto(data.pol())
    #     return core.snrRH(data.chop().data())

    def cmpangs(self):
        cmp1 = self.cmpvecs[:,0]
        cmp2 = self.cmpvecs[:,1]
        def getang(c) : return np.rad2deg(np.arctan2(c[1],c[0]))
        return getang(cmp1),getang(cmp2)          
    
    def chop(self):
        """
        Chop data to window
        """
        chop = self.copy()
        chop.x, chop.y = core.chop(chop.x,chop.y,window=chop.window)
        chop.window.offset = 0
        return chop
        
    def chopt(self):
        """
        Chop time to window
        """        
        t = core.chop(self.t(),window=self.window)
        return t
    
    def splitting_intensity(self,**kwargs):
        """
        Calculate the splitting intensity as defined by Chevrot (2000).
        """
        copy = self.copy()
        copy.rotateto(copy.pol)
        copy.x = np.gradient(copy.x)
        copy.chop()
        rdiff, trans = copy.x, copy.y
        s = -2 * np.trapz(trans * rdiff) / np.trapz(rdiff**2)
        return s
        
    # Window
    
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
        return (self.window.width-1) * self.delta

    def _construct_window(self,start,end,**kwargs): 
        if start > end: raise ValueError('start is larger than end')
        time_centre = (start + end)/2
        time_width = end - start
        tcs = core.time2samps(time_centre,self.delta)
        offset = tcs - self._centresamp()
        # convert time to nsamples -- must be odd (even plus 1 because x units of deltatime needs x+1 samples)
        width = core.time2samps(time_width,self.delta,'even') + 1     
        return Window(width,offset,**kwargs) 
        
    # Plotting
                
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
        ax1 = plt.subplot(gs[1])
        self._ppm( ax1, **kwargs)   
        
        # optional pick window
        if 'pick' in kwargs and kwargs['pick'] == True:
            windowpicker = WindowPicker(self,fig,ax0)
            windowpicker.connect()
                                 
        # show
        plt.tight_layout()
        plt.show()
        
    def ppm(self,**kwargs):
        """Plot particle motion"""
        fig, ax = plt.subplots()
        self._ppm(ax, **kwargs)
        plt.show()
        
    def ptr(self,**kwargs):
        """Plot trace data"""
        fig, ax = plt.subplots()
        self._ptr(ax, **kwargs)
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
        ax.legend(framealpha=0.5)
    
        # set limits
        lim = np.abs(self.data()).max() * 1.1
        if 'ylim' not in kwargs: kwargs['ylim'] = [-lim,lim]
        ax.set_ylim(kwargs['ylim'])
        if 'xlim' in kwargs: ax.set_xlim(kwargs['xlim'])
    
        # set axis label
        if 'units' not in kwargs: kwargs['units'] = 's'            
        ax.set_xlabel('Time (' + kwargs['units'] +')')

        # plot window markers
        if self.window.width < self._nsamps():
            w1 = ax.axvline(self.wbeg(),linewidth=1,color='k')
            w2 = ax.axvline(self.wend(),linewidth=1,color='k')    
        
        # plot additional markers
        if 'marker' in kwargs:
            print('here')
            if type(kwargs['marker']) is not list: kwargs['marker'] = [ kwargs['marker'] ]
            [ ax.axvline(float(mark),linewidth=1,color='b') for mark in kwargs['marker'] ]
            
        return

    def _ppm(self,ax,**kwargs):
        """Plot particle motion on *ax* matplotlib axis object.
        """
        
        data = self.chop()
        data.rotateto(0)
        x, y = data.x, data.y
        t = data.t()
                
        # plot data
        # ax.plot(self.chop().y,self.chop().x)
        
        # multi-colored
        norm = plt.Normalize(t.min(),t.max())
        points = np.array([y, x]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments,cmap='plasma',norm=norm,alpha=0.7)
        lc.set_array(t)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        # plt.colorbar(line)
    
        # set limit
        lim = np.abs(self.data()).max() * 1.1
        if 'lims' not in kwargs: kwargs['lims'] = [-lim,lim] 
        ax.set_aspect('equal')
        ax.set_xlim(kwargs['lims'])
        ax.set_ylim(kwargs['lims'])
    
        # set labels
        if 'cmplabels' not in kwargs: kwargs['cmplabels'] = data.cmplabels
        ax.set_xlabel(kwargs['cmplabels'][1])
        ax.set_ylabel(kwargs['cmplabels'][0])
        
        # turn off tick annotation
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
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

    return x,y
    


class WindowPicker:
    """
    Pick a Window
    """

    def __init__(self,pair,fig,ax):
           
        self.canvas = fig.canvas
        self.ax = ax
        self.pair = pair
        # window limit lines
        self.x1 = pair.wbeg()
        self.x2 = pair.wend()
        self.wbegline = self.ax.axvline(self.x1,linewidth=1,color='r',visible=True)
        self.wendline = self.ax.axvline(self.x2,linewidth=1,color='r',visible=True)
        self.cursorline = self.ax.axvline(pair._centretime(),linewidth=1,color='0.5',visible=False)
        _,self.ydat = self.wbegline.get_data()
            

    def connect(self):  
        self.cidclick = self.canvas.mpl_connect('button_press_event', self.click)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion)
        # self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        self.cidenter = self.canvas.mpl_connect('axes_enter_event', self.enter)
        self.cidleave = self.canvas.mpl_connect('axes_leave_event', self.leave)
        self.cidkey = self.canvas.mpl_connect('key_press_event', self.keypress) 
       
    def click(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        if event.button == 1:
            self.x1 = x
            self.wbegline.set_data([x,x],self.ydat)
            self.canvas.draw() 
        if event.button == 3:
            self.x2 = x
            self.wendline.set_data([x,x], self.ydat)
            self.canvas.draw()
    
    def keypress(self,event):
        if event.key == " ":
            self.disconnect()

    def enter(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.cursorline.set_visible(True)
        self.canvas.draw()

    def leave(self,event):
        if event.inaxes is not self.ax: return
        self.cursorline.set_visible(False)
        self.canvas.draw()

    def motion(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.canvas.draw()
        
    def disconnect(self):
        'disconnect all the stored connection ids'
        self.canvas.mpl_disconnect(self.cidclick)
        self.canvas.mpl_disconnect(self.cidmotion)
        self.canvas.mpl_disconnect(self.cidenter)
        self.canvas.mpl_disconnect(self.cidleave)
        plt.close()
        wbeg, wend = sorted((self.x1, self.x2)) 
        self.pair.set_window(wbeg, wend)
      