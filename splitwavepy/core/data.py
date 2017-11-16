# -*- coding: utf-8 -*-
"""
The seismic data class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core, core3d, io
# from ..core.pair import Pair
from ..core.window import Window
# from . import eigval, rotcorr, transmin, sintens

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os.path


class Data:
    
    """
    Base data class        
    """
    
    def __init__(self,*args,**kwargs):

        # ensure delta is set as a keyword argment, e.g. delta=0.1
        if 'delta' not in kwargs: raise Exception('delta must be set')
        self.delta = kwargs['delta']                    
        
    # COMMON PROPERTIES
    
    @property
    def delta(self):
        return self.__delta    

    @delta.setter
    def delta(self, delta):
        if delta <= 0: raise ValueError('delta must be positive')
        self.__delta = float(delta)
        
    @property
    def window(self):
        return self.__window
        
    @window.setter
    def window(self, window):    
        self.__window = window
   
    @property
    def units(self):
        return self.__units
    
    @units.setter
    def units(self, units):
        self.__units = units
                
    # Common methods
                
    def set_window(self,*args,**kwargs):
        """
        Set the window
        """                
        # if Window provided
        if 'window' in kwargs:  
            if isinstance(kwargs['window'],Window):
                self.window = kwargs['window']
            else:
                raise TypeError('expecting a window')        
        # if no arguments provided
        elif len(args) == 0:
            width = core.odd(self._nsamps() / 3)
            self.window = Window(width)
        # if start/end given
        elif len(args) == 2:
            start, end = args  
            self.window = self.construct_window(start,end,**kwargs) 
        else:
            raise Exception ('unexpected number of arguments')
            
        
    # Utility 
    
    def t(self):
        return np.arange(self._nsamps()) * self.delta
        
    def chopt(self):
        """
        Chop time to window
        """        
        t = core.chop(self.t(),window=self.window)
        return t
        
    # window
    
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
        
    def wcentre(self):
        """
        Window centre
        """
        return self.window.centre(self._nsamps()) * self.delta
        
    def construct_window(self,start,end,**kwargs): 
        if start > end: raise ValueError('start is larger than end')
        time_centre = (start + end)/2
        time_width = end - start
        tcs = core.time2samps(time_centre,self.delta)
        offset = tcs - self._centresamp()
        # convert time to nsamples -- must be odd (even plus 1 because x units of deltatime needs x+1 samples)
        width = core.time2samps(time_width,self.delta,'even') + 1     
        return Window(width,offset,**kwargs) 
        
    # Hidden
    
    def _nsamps(self):
        return self.x.size

    def _centresamp(self):
        return int(self.x.size/2)
    
    def _centretime(self):
        return int(self.x.size/2) * self.delta
                
    # I/O stuff  
                       
    def copy(self):
        return io.copy(self)

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
        
class WindowPicker:
    """
    Pick a Window
    """

    def __init__(self,data,fig,ax):
           
        self.canvas = fig.canvas
        self.ax = ax
        self.data = data
        # window limit lines
        self.x1 = data.wbeg()
        self.x2 = data.wend()
        self.wbegline = self.ax.axvline(self.x1,linewidth=1,color='r',visible=True)
        self.wendline = self.ax.axvline(self.x2,linewidth=1,color='r',visible=True)
        self.cursorline = self.ax.axvline(data._centretime(),linewidth=1,color='0.5',visible=False)
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
        self.data.set_window(wbeg, wend)