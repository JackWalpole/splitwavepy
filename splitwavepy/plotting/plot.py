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

import numpy as np

#### rcParams ####
# mpl.rcParams['axes.titlesize'] = 24
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['axes.titlepad'] = 12.0


def trace(*args,**kwargs):
    """Return axis with trace data.
    args:
    - x
    - y
    - z (optional)
    
    kwargs:
    - time
    - ax
    """
    
    if 'time' not in kwargs:
        kwargs['time'] = np.arange(args[0].size)   
    
    if 'ax' not in kwargs:        
        # initiate axis  
        kwargs['ax'] = plt.subplot(111)    
    ax = kwargs['ax']
        
    # plot data
    for ii in range(len(args)):
        ax.plot(kwargs['time'],args[ii])
        
    # set limit
    if 'ylim' not in kwargs:
        lim = np.abs(args).max() * 1.1
        kwargs['ylim'] = [-lim,lim]
    
    ax.set_ylim(kwargs['ylim'])    
        
    # set label
    if 'units' not in kwargs:
        kwargs['units'] = 's'  
          
    ax.set_xlabel('Time (' + kwargs['units'] +')')  
    
    # plot window markers
    if 'window' in kwargs:
        nsamps = args[0].size
        wbeg = kwargs['window'].start(nsamps)*kwargs['time'][1]
        wend = kwargs['window'].end(nsamps)*kwargs['time'][1]
        ax.axvline(wbeg,linewidth=1,color='k')
        ax.axvline(wend,linewidth=1,color='k')        
    
    return ax
    

def particle(*args,**kwargs):
    """Return axis with particle motion data
    
    kwargs:
    - ax
    """
    if not ('labels' in kwargs):
        kwargs['labels'] = ['comp1','comp2','comp3']
    
    # 2D particle motion
    if len(args) == 2:
        if 'ax' not in kwargs:         
            kwargs['ax'] = plt.subplot(111)
        ax = kwargs['ax']
        ax.plot(args[1],args[0])
        # set limit
        if 'lim' not in kwargs:
            lim = np.abs(args).max() * 1.1
            kwargs['lim'] = [-lim,lim]
        ax.set_aspect('equal')
        ax.set_xlim(kwargs['lim'])
        ax.set_ylim(kwargs['lim'])
        ax.set_xlabel(kwargs['labels'][1])
        ax.set_ylabel(kwargs['labels'][0])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        # ax.grid()
        return ax
    
    # 3D particle motion
    if len(args) == 3:
        if 'ax' not in kwargs:
            kwargs['ax'] = plt.subplot(111,projection='3d')
        ax = kwargs['ax']
        ax.plot(args[0],args[1],args[2])
        lim = np.abs(args).max() * 1.1
        ax.set_aspect('equal')
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])
        ax.set_zlim([-lim,lim])
        ax.plot(args[0],args[1],-lim,zdir='z',alpha=0.3,color='g')
        ax.plot(args[0],args[2],lim,zdir='y',alpha=0.3,color='g')
        ax.plot(args[1],args[2],-lim,zdir='x',alpha=0.3,color='g')
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
                
        # mpl.rcParams['keymap.back'].remove('c')
        # mpl.rcParams['keymap.all_axes'].remove('a')
        # mpl.rcParams['keymap.fullscreen'].remove('f')
        
        
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