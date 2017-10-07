# -*- coding: utf-8 -*-
"""
Some plotting routines
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# matplotlib stuff
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D

# splitwavepy stuff
# from ..core import core
from ..core import pair.Pair

# other stuff
import numpy as np

#### rcParams ####
# mpl.rcParams['axes.titlesize'] = 24
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['axes.titlepad'] = 12.0


def trace(*args,**kwargs):
    """Return axis with trace data.

    args:
    -Pair/Trio 
    
    kwargs:
    - ax = mpl axis
    """
    
    # Pairs or Trios only
    if not isinstance(args[0],Pair): 
        if not instance(args[0],Trio):
            raise(TypeError)
    
    # suck in the argument
    a = args[0]
        
    # accept an axis, or initiate or  
    if 'ax' not in kwargs: kwargs['ax'] = plt.subplot(111)          
    ax = kwargs['ax']
        
    # plot data
    ax.plot(a.t(),a.data())
        
    # set limit
    lim = np.abs(args).max() * 1.1
    if 'ylim' not in kwargs: kwargs['ylim'] = [-lim,lim]   
    ax.set_ylim(kwargs['ylim'])    
        
    # set label
    if 'units' not in kwargs: kwargs['units'] = 's'           
    ax.set_xlabel('Time (' + kwargs['units'] +')')  
    
    # plot window markers
    if 'window' in kwargs:
        ax.axvline(a.window.start(),linewidth=1,color='k')
        ax.axvline(a.window.end(),linewidth=1,color='k')        
    
    return ax
    

def particle(*args,**kwargs):
    """Return axis with particle motion data
    
    kwargs:
    - ax
    """
    
    # Pairs or Trios only
    if not isinstance(args[0],Pair): 
        if not instance(args[0],Trio):
            raise(TypeError)
    
    # suck in the argument
    a = args[0]
        

    
    if isinstance(a,Pair):
    
        # accept an axis, or initiate new 
        if 'ax' not in kwargs: kwargs['ax'] = plt.subplot(111)          
        ax = kwargs['ax']
    
        ax.plot(a.y,a.x)
    
        # set limit
        lim = np.abs(args).max() * 1.1
        if 'lim' not in kwargs: kwargs['lim'] = [-lim,lim]
    
        # axis properties
        ax.set_aspect('equal')
        ax.set_xlim(kwargs['lim'])
        ax.set_ylim(kwargs['lim'])
        ax.set_xlabel(kwargs['labels'][1])
        ax.set_ylabel(kwargs['labels'][0])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])

        return ax
    
    # 3D particle motion
    if isinstance(a,Trio):
    
        if 'ax' not in kwargs:
            kwargs['ax'] = plt.subplot(111,projection='3d')
        ax = kwargs['ax']
        
        # main data
        ax.plot(a.x,a.y,a.z)
        
        # fix limits
        lim = np.abs(a.all()).max() * 1.1
        ax.set_aspect('equal')
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
        ax.set_zlim([-lim,lim])
        
        # side panels
        ax.plot(a.x,a.y,-lim,zdir='z',alpha=0.3,color='g')
        ax.plot(a.x,a.z,lim,zdir='y',alpha=0.3,color='g')
        ax.plot(a.y,a.z,-lim,zdir='x',alpha=0.3,color='g')
        
        # axis properties
        ax.set_xlabel(kwargs['labels'][0])
        ax.set_ylabel(kwargs['labels'][1])
        ax.set_zlabel(kwargs['labels'][2])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.axes.zaxis.set_ticklabels([])
        
        return ax
        
        
def surf(M,**kwargs):
    """
    Plot an error surface.
    
    **kwargs
    - cmap = 'magma'
    - vals = (M.lam1-M.lam2) / M.lam2
    - ax = None (creates new)
    """
    
    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'magma'
    
    if 'vals' not in kwargs:
        kwargs['vals'] = (M.lam1-M.lam2) / M.lam2
            
    if 'ax' not in kwargs:         
        kwargs['ax'] = plt.subplot(111)
    
    ax = kwargs['ax']
        
    # error surface
    cax = ax.contourf(M.lags,M.degs,kwargs['vals'],26,cmap=kwargs['cmap'])
    cbar = plt.colorbar(cax)
    ax.set_yticks(np.linspace(-90,90,6,endpoint=False))
    ax.set_ylabel('Fast Direction (degs)')
    ax.set_xlabel('Delay Time (' + M.units + ')')
    
    # marker
    ax.errorbar(M.lag,M.fast,xerr=M.fdlag,yerr=M.fdfast,fmt='o')
    
    # confidence region
    ax.contour(M.lags,M.degs,M.lam2,levels=[M.lam2_95()])

    ax.set_xlim([M.lags[0,0], M.lags[-1,0]])
    ax.set_ylim([M.degs[0,0], M.degs[0,-1]])
    
    # optional title
    if 'title' in kwargs:
        ax.set_title(kwargs['title'])

    return ax


# Interaction
 
class WindowPicker:
    """
    Pick a Window
    """
    
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
        
        
        ### free up the keys I want to use so pyplot doesn't do funky things.
        ### probably a neater way to do this?
        neededkeys=['c','a','f',' ','enter']
        keymap = dict(plt.rcParams.find_all('keymap'))
        for key in keymap.keys():
            overlap = list(set(neededkeys) & set(keymap[key]))
            [ mpl.rcParams[key].remove(wantedkey) for wantedkey in overlap ]
        
        
        
        #
        # # the window
        # self.wbeg = None
        # self.wend = None
        # if 'window' in kwargs:
        #     self.wbeg = kwargs['window'].start()
        #     self.wend = kwargs['window'].end()

    def connect(self):
        # mouse interaction
        self.cidclick = self.canvas.mpl_connect('button_press_event', self.click)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion)
        self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        self.cidenter = self.canvas.mpl_connect('axes_enter_event', self.enter)
        self.cidleave = self.canvas.mpl_connect('axes_leave_event', self.leave)
        # # keyboard interaction
        self.cidkeypress = self.canvas.mpl_connect('key_press_event', self.keypress)
        # self.cidkeyrelease = self.canvas.mpl_connect('key_release_event', self.keyrelease)
    
               
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

    def keypress(self,event):
        
        if event.inaxes is not self.ax: return
        if event.key == 'c':
            self.wbegline.set_visible(False)
            self.wendline.set_visible(False)
        # establish which line is which
        xbeg = self.wbegline.get_data()[0][0]
        xend = self.wendline.get_data()[0][0]
        x = event.xdata
        if xbeg < xend: 
            leftline = self.wbegline
            rightline = self.wendline
        else:
            leftline = self.wendline
            rightline = self.wbegline
        if event.key == 'a':
            leftline.set_data([x,x],self.origydat)
            leftline.set_visible(True)
        if event.key == 'f':
            rightline.set_data([x,x],self.origydat)
            rightline.set_visible(True)
        if event.key == ' ':
            # chop the data and replot particle motion
            return 'halleluja!'
        # if event.key == 'q':
            # save window and quit

            
        
    # def keyrelease(self,event):
    #     if event.inaxes is self.ax: return
            
    # def disconnect(self):
    #     'disconnect all the stored connection ids'
    #     self.canvas.mpl_disconnect(self.cidclick)
    #     # self.canvas.mpl_disconnect(self.cidrelease)
    #     # self.canvas.mpl_disconnect(self.cidmotion)