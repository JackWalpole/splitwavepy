# -*- coding: utf-8 -*-
"""
The seismic data class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core, io

#, core3d, io
# from ..core.pair import Pair
# from ..core.window import Window
# from . import eigval, rotcorr, transmin, sintens

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection

class Data:
    
    """
    Base data class        
    """
    
    def __init__(self, x, y, *args, **kwargs):

        # the traces
        self.x = x
        self.y = y

        # ensure delta is set as a keyword argment, e.g. delta=0.1
        if 'delta' not in kwargs: raise Exception('delta must be set')
        self.delta = kwargs['delta'] 
        
        # some sanity checks
        if self.x.ndim != 1: raise Exception('data must be one dimensional')
        if (self.x.size != self.y.size): raise Exception('x and y must be the same length') 
        if self.x.size%2 == 0: 
            # drop last sample to ensure traces have odd number of samples
            self.x = x[:-1]
            self.y = y[:-1] 
        
        # add geometry info 
        self.geom = 'geo'
        if ('geom' in kwargs): self.geom = kwargs['geom']
               
        self.cmpvecs = np.eye(2)  
        if ('cmpvecs' in kwargs): self.cmpvecs = kwargs['cmpvecs']  
        
        # labels
        self.units = 's'
        if ('units' in kwargs): self.units = kwargs['units']
        self.set_labels()
        
        # default window
        self.window = Window(core.odd(self._nsamps() / 3))
                      
        
    # COMMON PROPERTIES
    
    
    
    @property
    def delta(self):
        return self.__delta    

    @delta.setter
    def delta(self, delta):
        if delta <= 0: raise ValueError('delta must be positive')
        self.__delta = float(delta)
        
    # @property
    # def window(self):
    #     return self.__window
    #
    # @window.setter
    # def window(self, window):
    #     self.__window = window
   
    @property
    def cmplabels(self):
        return self.__cmplabels
        
    @cmplabels.setter
    def cmplabels(self, cmplabels):
        self.__cmplabels = cmplabels
   
    @property
    def units(self):
        return self.__units
    
    @units.setter
    def units(self, units):
        self.__units = units
        
    @property
    def geom(self):
        return self.__geom
        
    @geom.setter
    def geom(self, geom):
        possible_geoms = ['geo', 'ray', 'cart']
        if geom not in possible_geoms:
            raise ValueError('geom must be one of ' + str(possible_geoms))
        self.__geom = geom
        
    # COMMON METHODS
    
    def split(self, fast, lag):
        """
        Applies splitting operator.
        
        .. warning:: shortens trace length by *lag*.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag, self.delta, mode='even')
        # find appropriate rotation angle
        origangs = self.cmpangs()
        self.rotateto(0)
        # apply splitting
        self.x, self.y = core.split(self.x, self.y, fast, samps)
        self.rotateto(origangs[0])
           
    def unsplit(self, fast, lag):
        """
        Reverses splitting operator.
        
        .. warning:: shortens trace length by *lag*.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag, self.delta, mode='even')
        # find appropriate rotation angle
        origangs=self.cmpangs()
        self.rotateto(0)
        # apply splitting
        self.x, self.y = core.unsplit(self.x, self.y, fast, samps)
        self.rotateto(origangs[0])
       
    def rotateto(self, degrees):
        """
        Rotate traces so that cmp1 lines up with *degrees*
        """
        # find appropriate rotation matrix
        ang = math.radians(degrees)
        cang = math.cos(ang)
        sang = math.sin(ang)
        # define the new cmpvecs
        backoff = self.cmpvecs
        self.cmpvecs = np.array([[ cang,-sang],
                                 [ sang, cang]])
        rot = np.dot(self.cmpvecs.T, backoff)
        # rotate data
        xy = np.dot(rot, self.data())
        self.x, self.y = xy[0], xy[1]
        # reset label
        self.set_labels()
                
    def set_window(self, *args, **kwargs):
        """
        Set the window
        """                
        # if Window provided
        if 'window' in kwargs:  
            if isinstance(kwargs['window'], Window):
                self.window = kwargs['window']
                return
            else:
                raise TypeError('expecting a window')  
        # start/end given
        if len(args) == 2:
            start, end = args  
            self.window = self.construct_window(start, end, **kwargs)
            return
        else:
            raise Exception ('unexpected number of arguments')
            
    def set_labels(self, *args):
        if len(args) == 0:
            if np.allclose(self.cmpvecs, np.eye(2), atol=1e-02):
                if self.geom == 'geo': self.cmplabels = ['North', 'East']
                elif self.geom == 'ray': self.cmplabels = ['SV', 'SH']
                elif self.geom == 'cart': self.cmplabels = ['X', 'Y']
                else: self.cmplabels = ['Comp1', 'Comp2']
                return
            # if reached here we have a non-standard orientation
            a1,a2 = self.cmpangs()
            lbl1 = str(round(a1))+r' ($^\circ$)'
            lbl2 = str(round(a2))+r' ($^\circ$)'
            self.cmplabels = [lbl1, lbl2]
            return
        elif len(args) == 1:
            if not isinstance(args[0], list): raise TypeError('expecting a list')
            # if not len(args[0]) == 2: raise Exception('list must be length 2')
            if not (isinstance(args[0][0], str) and isinstance(args[0][1], str)):
                raise TypeError('cmplabels must be a list of strings')
            self.cmplabels = args[0]
            return
        else:
            raise Exception('unexpected number of arguments')     
        
    # Utility 
    
    def t(self):
        return np.arange(self._nsamps()) * self.delta
        
    def chopt(self):
        """
        Chop time to window
        """
        t0 = self._w0()
        t1 = self._w1()        
        t = self.t()[t0:t1]
        return t
        
    def data(self):
        return np.vstack((self.x, self.y))
        
    def chopdata(self):
        """Chop traces to window"""
        t0 = self._w0()
        t1 = self._w1()
        return self.x[t0:t1], self.y[t0:t1]
        # return np.vstack((self.x[t0:t1], self.y[t0:t1]))
        
    def chop(self):
        chop = self.copy()
        chop.x, chop.y = chop.chopdata()
        chop.window.offset = 0
        return chop
        
    def estimate_pol(self):
        """Return principal component orientation"""
        # rotate to zero
        rot = self.cmpvecs.T
        data = np.vstack((self.chopdata()))
        xy = np.dot(rot,data)
        _,eigvecs = core.eigcov(xy)
        x,y = eigvecs[:,0]
        pol = np.rad2deg(np.arctan2(y,x))
        return pol
        
    # window
    
    def _w0(self):
        """idx of first sample in window"""
        hw = int(self.window.width/2)
        return self._centresamp() + self.window.offset - hw
    
    def _w1(self):
        """idx of last sample in window"""
        hw = int(self.window.width/2)
        return self._centresamp() + self.window.offset + hw + 1
    
    def wbeg(self):
        """
        Window start time.
        """
        sbeg = self._w0()
        return sbeg * self.delta
    
    def wend(self):
        """
        Window end time.
        """
        send = self._w1()
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
        
    def construct_window(self, start, end, **kwargs): 
        if start > end: raise ValueError('start is larger than end')
        time_centre = (start + end)/2
        time_width = end - start
        tcs = core.time2samps(time_centre, self.delta)
        offset = tcs - self._centresamp()
        # convert time to nsamples -- must be odd (even plus 1 because x units of deltatime needs x+1 samples)
        width = core.time2samps(time_width, self.delta, 'even') + 1     
        return Window(width, offset, **kwargs) 
        
    def eigen(self, window=None):
        self.eigvals, self.eigvecs = core.eigcov(self.data())
        
    def power(self):
        return self.x**2, self.y**2
        
    # def snrRH(self):
    #     data = self.copy()
    #     data.rotateto(data.pol())
    #     return core.snrRH(data.chop().data())

    def cmpangs(self):
        cmp1 = self.cmpvecs[:, 0]
        cmp2 = self.cmpvecs[:, 1]
        def getang(c) : return np.rad2deg(np.arctan2(c[1], c[0]))
        return getang(cmp1), getang(cmp2)
        
    # Hidden
    
    
    # Plotting
              
    def plot(self, **kwargs):
        """
        Plot trace data and particle motion
        """

        fig = plt.figure(figsize=(12, 3))     
        gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
        
        # trace
        ax0 = plt.subplot(gs[0])
        self._ptr(ax0, **kwargs)
        
        # particle  motion
        ax1 = plt.subplot(gs[1])
        self._ppm( ax1, **kwargs)   
        
        # optional pick window
        if 'pick' in kwargs and kwargs['pick'] == True:
            windowpicker = WindowPicker(self, fig, ax0)
            windowpicker.connect()
                                 
        # show
        plt.tight_layout()
        plt.show()
        
    def ppm(self, **kwargs):
        """Plot particle motion"""
        fig, ax = plt.subplots()
        self._ppm(ax, **kwargs)
        plt.show()
        
    def ptr(self, **kwargs):
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
        if 'ylim' not in kwargs: kwargs['ylim'] = [-lim, lim]
        ax.set_ylim(kwargs['ylim'])
        if 'xlim' in kwargs: ax.set_xlim(kwargs['xlim'])
    
        # set axis label
        if 'units' not in kwargs: kwargs['units'] = 's'            
        ax.set_xlabel('Time (' + kwargs['units'] +')')

        # plot window markers
        if self.window.width < self._nsamps():
            w1 = ax.axvline(self.wbeg(), linewidth=1, color='k')
            w2 = ax.axvline(self.wend(), linewidth=1, color='k')    
        
        # plot additional markers
        if 'marker' in kwargs:
            print('here')
            if type(kwargs['marker']) is not list: kwargs['marker'] = [ kwargs['marker'] ]
            [ ax.axvline(float(mark), linewidth=1, color='b') for mark in kwargs['marker'] ]
            
        return

    def _ppm(self, ax, **kwargs):
        """Plot particle motion on *ax* matplotlib axis object.
        """
        
        data = self.copy()
        data.rotateto(0)
        x, y = data.chopdata()
        t = data.chopt()
                
        # plot data
        # ax.plot(self.chop().y,self.chop().x)
        
        # multi-colored
        norm = plt.Normalize(t.min(), t.max())
        points = np.array([y, x]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap='plasma', norm=norm, alpha=0.7)
        lc.set_array(t)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        # plt.colorbar(line)
    
        # set limit
        lim = np.abs(self.data()).max() * 1.1
        if 'lims' not in kwargs: kwargs['lims'] = [-lim, lim] 
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
        
class Window:
    """
    Instantiate a Window defined relative to centre of a window of flexible size.
    
    args

    - width    | nsamps length of window,
    - offset   | nsamps offset from centre of window,    
    
    kwargs
    
    - tukey   | fraction of window to cosine taper (from 0 to 1).
    """
    
    def __init__(self, width, offset=0, tukey=None):
        # ensure width is odd 
        if width%2 != 1:
            raise Exception('width must be an odd integer')
        self.width = width
        self.offset = offset
        self.tukey = tukey
    
    # def start(self, samps):
    #     """
    #     Return start sample of window.
    #     """
    #     hw = int(self.width/2)
    #     if samps%2 != 1:
    #         raise Exception('samps must be odd to have definite centre')
    #     else:
    #         centre = np.int(samps/2)
    #         return centre + self.offset - hw
    #
    # def end(self, samps):
    #     """
    #     Return end sample of window.
    #     """
    #     hw = int(self.width/2)
    #     if samps%2 != 1:
    #         raise Exception('samps must be odd to have definite centre')
    #     else:
    #         centre = int(samps/2)
    #         return centre + self.offset + hw
    #
    # def centre(self, samps):
    #     """
    #     Return centre sample of window.
    #     """
    #     if samps%2 != 1:
    #         raise Exception('samps must be odd to have definite centre')
    #     else:
    #         centre = int(samps/2)
    #         return centre + self.offset
    #
    # def asarray(self, samps):
    #
    #     # sense check -- is window in range?
    #     if self.end(samps) > samps:
    #         raise Exception('Window exceeds max range')
    #     if self.start(samps) < 0:
    #         raise Exception('Window exceeds min range')
    #
    #     # sexy cosine taper
    #     if self.tukey is None:
    #         alpha = 0.
    #     else:
    #         alpha = self.tukey
    #     tukey = signal.tukey(self.width, alpha=alpha)
    #     array = np.zeros(samps)
    #     array[self.start(samps):self.end(samps)+1] = tukey
    #     return array
                
    # def shift(self, shift):
    #     """
    #     +ve moves N samples to the right
    #     """
    #     self.offset = self.offset + int(shift)
    #
    # def resize(self, resize):
    #     """
    #     +ve adds N samples to the window width
    #     """
    #     # ensure resize is even
    #     self.width = self.width + core.even(resize)
    #
    # def retukey(self, tukey):
    #     self.tukey = tukey
        
    # Comparison
    
    def __eq__(self, other) :
        if self.__class__ != other.__class__: return False
        if set(self.__dict__) != set(other.__dict__): return False
        return True



class WindowPicker:
    """
    Pick a Window
    """

    def __init__(self, data, fig, ax):
           
        self.canvas = fig.canvas
        self.ax = ax
        self.data = data
        # window limit lines
        self.x1 = data.wbeg()
        self.x2 = data.wend()
        self.wbegline = self.ax.axvline(self.x1, linewidth=1, color='r', visible=True)
        self.wendline = self.ax.axvline(self.x2, linewidth=1, color='r', visible=True)
        self.cursorline = self.ax.axvline(data._centretime(), linewidth=1, color='0.5', visible=False)
        _, self.ydat = self.wbegline.get_data()
            
    def connect(self):  
        self.cidclick = self.canvas.mpl_connect('button_press_event', self.click)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion)
        # self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        self.cidenter = self.canvas.mpl_connect('axes_enter_event', self.enter)
        self.cidleave = self.canvas.mpl_connect('axes_leave_event', self.leave)
        self.cidkey = self.canvas.mpl_connect('key_press_event', self.keypress) 
       
    def click(self, event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        if event.button == 1:
            self.x1 = x
            self.wbegline.set_data([x, x], self.ydat)
            self.canvas.draw() 
        if event.button == 3:
            self.x2 = x
            self.wendline.set_data([x, x], self.ydat)
            self.canvas.draw()
    
    def keypress(self, event):
        if event.key == " ":
            self.disconnect()

    def enter(self, event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x, x], self.ydat)
        self.cursorline.set_visible(True)
        self.canvas.draw()

    def leave(self, event):
        if event.inaxes is not self.ax: return
        self.cursorline.set_visible(False)
        self.canvas.draw()

    def motion(self, event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x, x], self.ydat)
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