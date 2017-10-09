# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core
# from ..plotting import plot
from .window import Window
from ..core import io
from . import geom

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
        if ('cmpvecs' in kwargs): self.cmpvecs = kwargs['cmpvecs']
        
        self.rayvec = [0,0,1] # normal to shear plane, along Z-axis
        if ('rayvec' in kwargs): self.rayvec = kwargs['rayvec']
        
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

    # methods

    def t(self):
        return np.arange(self.x.size) * self.delta
  
    def data(self):
        return np.vstack((self.x,self.y))
    
    def set_labels(self):
        if np.allclose(self.cmpvecs,np.eye(2),atol=1e-02):
            if self.geom == 'geo': self.cmplabels = ['North','East']
            elif self.geom == 'ray': self.cmplabels = ['Vertical','Horizontal']
            elif self.geom == 'cart': self.cmplabels = ['X','Y']
            else: self.cmplabel = ['Comp1','Comp2']
            return
        # if reached here we have a non-standard orientation
        a1,a2 = self.cmpangs()
        lbl1 = str(round(a1))+r' ($^\circ$)'
        lbl2 = str(round(a2))+r' ($^\circ$)'
        self.cmplabels = [lbl1,lbl2]
        
        
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
        # reset label
        # if reached here use the default label
        lab1,lab2 = self.cmpangs()
        self.set_labels()
        
        
    def pol(self):
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

    def cmpangs(self):
        cmp1 = self.cmpvecs[:,0]
        cmp2 = self.cmpvecs[:,1]
        def getang(c) : return np.rad2deg(np.arctan2(c[1],c[0]))
        return getang(cmp1),getang(cmp2)
        
# Geometry stuff

    # def geom_to_geo():
    # def geom_to_ray():
    # def geom_to

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
            tcs = core.time2samps(time_centre,self.delta)
            offset = tcs - self.centresamp()
            # convert time to nsamples -- must be odd
            width = core.time2samps(time_width,self.delta,'odd')       
            self.window = Window(width,offset,**kwargs) 
            return
        else:
            raise Exception ('unexpected number of arguments')     
    
    def chop(self):
        """
        Chop data to window
        """
        chop = self.copy()
        chop.x, chop.y = core.chop(self.x,self.y,window=self.window)
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
        sbeg = self.window.start(self.nsamps())
        return sbeg * self.delta
    
    def wend(self):
        """
        Window end time.
        """
        send = self.window.end(self.nsamps())
        return send * self.delta
    
    # plotting
                
    def plot(self,**kwargs):
        """
        Plot trace data and particle motion
        """
        from matplotlib import gridspec
        fig = plt.figure(figsize=(12, 3))
        
        gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
        
        # trace
        ax0 = plt.subplot(gs[0])
        self._ptr( ax0, **kwargs)
        
        # particle  motion
        ax1 = plt.subplot(gs[1])
        self._ppm( ax1, **kwargs)
        
        if 'interact' in kwargs:
            windowpicker = WindowPicker(fig,ax0)
            windowpicker.connect()
                                 
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
        if self.window.width < self.nsamps():
            ax.axvline(self.wbeg(),linewidth=1,color='k',picker=100)
            ax.axvline(self.wend(),linewidth=1,color='k',picker=100)       

        return

    def _ppm(self,ax,**kwargs):
        """Plot particle motion on *ax* matplotlib axis object.
        """
                
        # plot data
        ax.plot(self.chop().y,self.chop().x)
    
        # set limit
        lim = np.abs(self.data()).max() * 1.1
        if 'lims' not in kwargs: kwargs['lims'] = [-lim,lim] 
        ax.set_aspect('equal')
        ax.set_xlim(kwargs['lims'])
        ax.set_ylim(kwargs['lims'])
    
        # set labels
        if 'cmplabels' not in kwargs: kwargs['cmplabels'] = self.cmplabels
        ax.set_xlabel(kwargs['cmplabels'][1])
        ax.set_ylabel(kwargs['cmplabels'][0])
        
        # turn off tick annotation
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])

        return ax
        
    # useful?
    def nsamps(self):
        return self.x.size
    
    def centresamp(self):
        return int(self.x.size/2)
        
    def centretime(self):
        return int(self.x.size/2) * self.delta
        
    def power(self):
        return self.x**2+self.y**2
      
    # I/O stuff  

    def save(self,filename):
        """
        Save Measurement for future referral
        """
        io.save(self,filename)
                       
    def copy(self):
        return copy.deepcopy(self)    

    # In builts
    
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
        
# Interaction

class WindowPicker:
    """
    Pick a Window
    """

    ### free up the keys I want to use so pyplot doesn't do funky things.
    ### probably a neater way to do this?
    # import matplotlib as mpl
    # neededkeys=['c','a','f',' ','enter']
    # keymap = dict(plt.rcParams.find_all('keymap'))
    # for key in keymap.keys():
    #     overlap = list(set(neededkeys) & set(keymap[key]))
    #     [ mpl.rcParams[key].remove(wantedkey) for wantedkey in overlap ]

    def __init__(self,fig,ax,**kwargs):

        self.canvas = fig.canvas
        self.ax = ax

        self.press = None

        # window limit lines
        minx, maxx = self.ax.get_xlim()
        self.wbegline = self.ax.axvline(minx,linewidth=1,color='k',visible=False)
        self.wendline = self.ax.axvline(maxx,linewidth=1,color='k',visible=False)
        self.cursorline = self.ax.axvline((minx+maxx)/2,linewidth=1,color='0.5',visible=False)
        _,self.origydat = self.wbegline.get_data()

    def onpick(self,event):
        ind = event.ind
        xdata = event.artist.get_xdata()
        print('this line x:', xdata,)
        msx = event.mouseevent.xdata
        dist = abs(msx-xdata)
        ind = [ind[np.argmin(dist)]]
        print(dist,ind)

            

    def connect(self):
        #
        # # mouse interaction
        # self.cidclick = self.canvas.mpl_connect('button_press_event', self.click)
        # self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion)
        # self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        # self.cidenter = self.canvas.mpl_connect('axes_enter_event', self.enter)
        # self.cidleave = self.canvas.mpl_connect('axes_leave_event', self.leave)
        # # # keyboard interaction
        # self.cidkeypress = self.canvas.mpl_connect('key_press_event', self.keypress)
        # # self.cidkeyrelease = self.canvas.mpl_connect('key_release_event', self.keyrelease)
        self.canvas.mpl_connect('pick_event',self.onpick)
       
    def click(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.wbegline.set_data([x,x],self.origydat)
        self.wbegline.set_visible(True)
        self.canvas.draw()

    def release(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.wendline.set_data([x,x],self.origydat)
        self.wendline.set_visible(True)
        self.canvas.draw()
        
    def enter(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.origydat)
        self.cursorline.set_visible(True)
        self.canvas.draw()

    def leave(self,event):
        if event.inaxes is not self.ax: return
        self.cursorline.set_visible(False)
        self.canvas.draw()

    def motion(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.origydat)
        self.canvas.draw()

    # def keypress(self,event):
    #
    #     if event.inaxes is not self.ax: return
    #     if event.key == 'c':
    #         self.wbegline.set_visible(False)
    #         self.wendline.set_visible(False)
    #     # establish which line is which
    #     xbeg = self.wbegline.get_data()[0][0]
    #     xend = self.wendline.get_data()[0][0]
    #     x = event.xdata
    #     if xbeg < xend:
    #         leftline = self.wbegline
    #         rightline = self.wendline
    #     else:
    #         leftline = self.wendline
    #         rightline = self.wbegline
    #     if event.key == 'a':
    #         leftline.set_data([x,x],self.origydat)
    #         leftline.set_visible(True)
    #     if event.key == 'f':
    #         rightline.set_data([x,x],self.origydat)
    #         rightline.set_visible(True)
    #     if event.key == ' ':
    #         # chop the data and replot particle motion
    #         return 'halleluja!'
    
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
    # core.rotate does the passive transormation of the co-ordinate system
    x,y = core.rotate(x,y,-pol)    
    # add any splitting -- lag samples must be even
    slag = core.time2samps(lag,delta,mode='even')
    x,y = core.split(x,y,fast,slag)
    
    return x,y
    


