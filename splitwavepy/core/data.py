# -*- coding: utf-8 -*-
"""
The seismic data class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
# from ..core.pair import Pair
# from ..core.window import Window
# from . import eigval, rotcorr, transmin, sintens

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os.path
from scipy import signal


class Data:
    
    """
    Base data class        
    """
    
    def __init__(self,*args,**kwargs):

        # ensure delta is set as a keyword argument, e.g. delta=0.1
        if 'delta' not in kwargs: raise Exception('delta must be set')
        self.delta = kwargs['delta']
        
        if len(args) == 2:
            self.x, self.y = args[0], args[1]
            self.__ncomps = 2
        elif len(args) == 3:
            self.x, self.y, self.z = args[0], args[1], args[2]
            self.__ncomps = 3
        else:
            raise Exception('expects x, y, (and maybe z) in numpy arrays')
            
        # defaults
        self.cmpvecs = np.eye(3)
        self.set_window()      
        
    # COMMON PROPERTIES
    
    @property
    def delta(self):
        return self.__delta    

    @delta.setter
    def delta(self, delta):
        if delta <= 0: raise ValueError('delta must be positive')
        self.__delta = float(delta)
        
    @property
    def pol(self):
        return self.__pol
        
    @pol.setter
    def pol(self, pol):
        self.__pol = float(pol)
        
    # @property
    # def polvec(self):
    #     return self.__polvec
    #
    # @polvec.setter
    # def polvec(self, polvec):
    #     if not isinstance(polvec, np.ndarray):
    #         raise TypeError('polvec must be a numpy array')
        
    @property
    def window(self):
        return self.__window
        
    @window.setter
    def window(self, window):    
        self.__window = window
        
    @property
    def cmpvecs(self):
        return self.__cmpvecs
        
    @cmpvecs.setter
    def cmpvecs(self, cmpvecs):
        if not isinstance(cmpvecs, np.ndarray):
            raise TypeError('cmpvecs must be a numpy array')
        # check components are orthogonal
        self.__cmpvecs = cmpvecs
   
    # @property
    # def rayvec(self):
    #     return self.__rayvec
    #
    # @rayvec.setter
    # def rayvec(self, rayvec):
    #     self.__rayvec = rayvec
        
   
    # @property
    # def units(self):
    #     return self.__units
    #
    # @units.setter
    # def units(self, units):
    #     self.__units = units
                
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
            
    def set_labels(self,*args):
        if len(args) == 0:
            if np.allclose(self.cmpvecs,np.eye(3),atol=1e-02):
                self.cmplabels = ['North','East','Up'] #
                return
            if np.allclose(self.rayvecs,self.cmpvecs):
                self.cmplabels = ['SV','SH','P']
                return
            # if reached here we have a non-standard orientation
            self.cmplabels = ['Comp1','Comp2','Comp3']
            return
        elif len(args) == 1:
            if not isinstance(args[0],list): raise TypeError('expecting a list')
            # if not len(args[0]) == 3: raise Exception('list must be length 3')
            if not (isinstance(args[0][0],str) and 
                    isinstance(args[0][1],str) and
                    isinstance(args[0][2],str)):
                raise TypeError('cmplabels must be a list of strings')
            self.cmplabels = args[0]
            return
        else:
            raise Exception('unexpected number of arguments')
    
    # def set_pol(self, **kwargs):
    #     """
    #     Set polarisation
    #     """
    #     if 'pol' not in kwargs:
    #         raise Exception('pol must be a keyword argument')
    #
    #     self.pol = kwargs['pol']
        
    # Time
    
    def t(self):
        return np.arange(self._nsamps()) * self.__delta
        
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
        return sbeg * self.__delta

    def wend(self):
        """
        Window end time.
        """
        send = self.window.end(self._nsamps())
        return send * self.__delta
        
    def wwidth(self):
        """
        Window width.
        """
        return (self.window.width-1) * self.__delta
        
    def wcentre(self):
        """
        Window centre
        """
        return self.window.centre(self._nsamps()) * self.__delta
        
    def construct_window(self,start,end,**kwargs): 
        if start > end: raise ValueError('start is larger than end')
        time_centre = (start + end)/2
        time_width = end - start
        tcs = core.time2samps(time_centre,self.__delta)
        offset = tcs - self._centresamp()
        # convert time to nsamples -- must be odd (even plus 1 because x units of deltatime needs x+1 samples)
        width = core.time2samps(time_width,self.__delta,'even') + 1     
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
        
# import numpy as np
# import matplotlib.pyplot as plt

class Window:
    """
    Instantiate a Window defined relative to centre of a window of flexible size.
    
    args

    - width    | nsamps length of window,
    - offset   | nsamps offset from centre of window,    
    
    kwargs
    
    - tukey   | fraction of window to cosine taper (from 0 to 1).
    """
    
    def __init__(self,width,offset=0,tukey=None):
        # ensure width is odd 
        if width%2 != 1:
            raise Exception('width must be an odd integer')
        self.width = width
        self.offset = offset
        self.tukey = tukey
    
    def start(self,samps):
        """
        Return start sample of window.
        """
        hw = int(self.width/2)
        if samps%2 != 1:
            raise Exception('samps must be odd to have definite centre')
        else:
            centre = np.int(samps/2)
            return centre + self.offset - hw

    def end(self,samps):
        """
        Return end sample of window.
        """
        hw = int(self.width/2)
        if samps%2 != 1:
            raise Exception('samps must be odd to have definite centre')
        else:
            centre = int(samps/2)
            return centre + self.offset + hw
    
    def centre(self,samps):
        """
        Return centre sample of window.
        """
        if samps%2 != 1:
            raise Exception('samps must be odd to have definite centre')
        else:
            centre = int(samps/2)
            return centre + self.offset       

    def asarray(self,samps):
                
        # sense check -- is window in range?
        if self.end(samps) > samps:
            raise Exception('Window exceeds max range')        
        if self.start(samps) < 0:
            raise Exception('Window exceeds min range')
        
        # sexy cosine taper
        if self.tukey is None:
            alpha = 0.
        else:
            alpha = self.tukey
        tukey = signal.tukey(self.width,alpha=alpha)        
        array = np.zeros(samps)
        array[self.start(samps):self.end(samps)+1] = tukey
        return array
                
    def shift(self,shift):
        """
        +ve moves N samples to the right
        """
        self.offset = self.offset + int(shift)
        
    def resize(self,resize):
        """
        +ve adds N samples to the window width
        """        
        # ensure resize is even
        self.width = self.width + core.even(resize)
        
    def retukey(self,tukey):
        self.tukey = tukey
        
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