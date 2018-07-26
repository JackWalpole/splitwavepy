#!/usr/bin/env python

import splitwavepy as sw
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection
import numpy as np

# def plot(self, **kwargs):
#     """
#     Plot trace data and particle motion
#     """
#
#     settings = {}
#     settings['pick'] = False
#     settings.update(**kwargs)
#
#     fig = plt.figure(figsize=(12, 3), **kwargs)
#     gs = gridspec.GridSpec(1, 3, width_ratios=[3, 1, 2], **kwargs)
#
#     # trace
#     ax0 = plt.subplot(gs[0])
#     linex, liney = _ptr(self._data, ax0, **kwargs)
#
#     # particle  motion
#     ax1 = plt.subplot(gs[1])
#     _ppm(self._data, ax1, **kwargs)
#
#     # surface
#     ax2 = plt.subplot(gs[2])
#     _psurf(self, ax2, vals=self.sc.vals, **kwargs)
#     #
#     # # optional pick window
#     # if settings['pick'] == True:
#     #     windowpicker = WindowPicker(self, fig, ax0)
#     #     windowpicker.connect()
#
#     # explore
#     explorer = Explorer(self, fig, ax0, ax1, ax2)
#     explorer.connect()
#
#
#     # neaten
#     plt.tight_layout()
#
#     plt.show()
        

        
def _ptr(self, ax, **kwargs):
    """Plot trace data on *ax* matplotlib axis object.
    """

    # plot data
    t = self._t
    linex = ax.plot( t, self.x, label=self._labels[0])
    liney = ax.plot( t, self.y, label=self._labels[1])
    lim = np.abs(self._xy).max() * 1.1
    ax.set_ylim([-lim, lim])
    ax.set_xlabel('Time (' + self.units +')')
    ax.legend(framealpha=0.5)

    # plot window markers
    if self._window.width < self._nsamps():
        w1 = ax.axvline(self.wbeg(), linewidth=1, color='k')
        w2 = ax.axvline(self.wend(), linewidth=1, color='k')    
    
    # plot additional markers
    if 'marker' in kwargs:
        if type(kwargs['marker']) is not list: kwargs['marker'] = [ kwargs['marker'] ]
        [ ax.axvline(float(mark), linewidth=1, color='b') for mark in kwargs['marker'] ]
        
    return linex[0], liney[0], w1, w2

def _ppm(self, ax, **kwargs):
    """Plot particle motion on *ax* matplotlib axis object.
    """
    
    data = self.rotateto(0)
    x, y = data._chopxy()
    t = data._chopt()
            
    # plot data
    lines = ax.plot(y, x)
    
    # multi-colored
    # norm = plt.Normalize(t.min(), t.max())
    # points = np.array([y, x]).T.reshape(-1, 1, 2)
    # segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # lc = LineCollection(segments, cmap='plasma', norm=norm, alpha=0.7)
    # lc.set_array(t)
    # lc.set_linewidth(2)
    # line = ax.add_collection(lc)
    # plt.colorbar(line)

    # set limit
    lim = np.abs(self._xy).max() * 1.1
    if 'lims' not in kwargs: kwargs['lims'] = [-lim, lim] 
    ax.set_aspect('equal')
    ax.set_xlim(kwargs['lims'])
    ax.set_ylim(kwargs['lims'])

    # set labels
    if 'label' not in kwargs: kwargs['label'] = data._labels
    ax.set_xlabel(kwargs['label'][1])
    ax.set_ylabel(kwargs['label'][0])
    
    # turn off tick annotation
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    
    return lines[0]

def _psurf(self, ax, **kwargs):
    """
    Plot an error surface.

    **kwargs
    - cmap = 'magma'
    - vals = (M.lam1-M.lam2) / M.lam2
    """

    if 'cmap' not in kwargs:
        kwargs['cmap'] = 'magma'

    if 'vals' not in kwargs:
        raise Exception('vals must be specified')
        
    laggrid, deggrid = self._grid
        
    # error surface
    cax = ax.contourf(laggrid, deggrid, kwargs['vals'], 26, cmap=kwargs['cmap'])
    # cbar = plt.colorbar(cax)
    ax.set_ylabel(r'Fast Direction ($^\circ$)')
    ax.set_xlabel('Delay Time (' + self._data.units + ')')
    
    # confidence region
    # if 'conf95' in kwargs and kwargs['conf95'] == True:
    #     ax.contour(self.lagmap, self.degmap, self.errsurf, levels=[self.conf95level],
    #             colors='r', alpha=.5, linewidths=3)
        
    # marker
    if 'marker' in kwargs and kwargs['marker'] == True:
        ax.errorbar(self.lag, self.fast, xerr=self.dlag, yerr=self.dfast)

    ax.set_xlim([laggrid[0,0], laggrid[-1,0]])
    ax.set_ylim([deggrid[0,0], deggrid[0,-1]])

    # optional title
    if 'title' in kwargs:
        ax.set_title(kwargs['title']) 
        
    # add info in text box
    if 'info' in kwargs and kwargs['info'] == True:
        textstr = '$\phi=%.1f\pm%.1f$\n$\delta t=%.2f\pm%.2f$'%\
                    (self.fast, self.dfast, self.lag, self.dlag)
        # place a text box in upper left in axes coords
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=props)
    
    if 'ppm' in kwargs and kwargs['ppm'] is True:
        sublags = self.lags[0:-1:int(self.lags.size/6)]
        subdegs = self.degs[0:-1:int(self.degs.size/6)]
        sublags = sublags + (self.lags[-1]-sublags[-1]) / 2
        subdegs = subdegs + (self.degs[-1]-subdegs[-1]) / 2
        x, y = self.SplitWave_corr()._chopdata()   
        lagtot = self.lags[-1] - self.lags[0]
        degtot = self.degs[-1] - self.degs[0]
        boost = 10 * lagtot / np.max((x**2 + y**2)**.5)      
        for fast in subdegs:
            for lag in sublags:
                x, y = self.unsplit(fast, lag)._chopdata()
                ax.plot(lag + y*boost/degtot, fast + x*boost/lagtot, color='w',alpha=0.5)
    
    return cax
    

class Explorer:
    
    def __init__(self, Py, **kwargs):
        # self.m = m
        # self.fig = fig
        # self.ax0 = ax0
        # self.ax1 = ax1
        # self.ax2 = ax2
        
        self.py = Py

        self.fig = plt.figure(figsize=(18, 3), **kwargs)     
        gs = gridspec.GridSpec(1, 4, width_ratios=[5, 1, 2, 2], **kwargs) 
        
    
        # trace
        self.ax0 = plt.subplot(gs[0])
        self.linex, self.liney, self.w1, self.w2 = _ptr(Py._data, self.ax0, **kwargs)
        _, self.ydat = self.w1.get_data()
    
        # particle  motion
        self.ax1 = plt.subplot(gs[1])
        self.ppm = _ppm(Py._data, self.ax1, **kwargs)
    
        # surface silver and chan
        self.ax2 = plt.subplot(gs[2])
        self.sc_surf = _psurf(self.py, self.ax2, vals=Py.sc.vals, **kwargs)  
        
        # surface cross-correlation
        self.ax3 = plt.subplot(gs[3])
        self.xc_surf = _psurf(self.py, self.ax3, vals=Py.xc.vals, **kwargs) 
         
        #
        # # optional pick window
        # if settings['pick'] == True:
        #     windowpicker = WindowPicker(self, fig, ax0)
        #     windowpicker.connect()

        # explore
        # explorer = Explorer(self, fig, ax0, ax1, ax2)
        self.connect()

                         
        # neaten
        plt.tight_layout()

        plt.show()
        
    def connect(self):
        self.on_hover = self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.on_click = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def on_hover(self, event):
        """hover over surface"""
        lag, fast = event.xdata, event.ydata
        if (event.inaxes == self.ax2) | (event.inaxes == self.ax3):
            d = self.py._data.unsplit(fast, lag)
            self.linex.set_data(d._t, d.x)
            self.liney.set_data(d._t, d.y)
            chopxy = d._chopxy()
            self.ppm.set_data(chopxy[1], chopxy[0])
        plt.draw()        
            
    def on_click(self, event, **kwargs):
        """click trace window"""
        x = event.xdata
        if event.inaxes == self.ax0:
            # left click
            if event.button == 1:
                self.w1.set_data([x, x], self.ydat)
            # right click
            elif event.button == 3:
                self.w2.set_data([x, x], self.ydat)
            else: return
            
            x1, _ = self.w1.get_data()
            x2, _ = self.w2.get_data()
            wbeg, wend = sorted((x1[0], x2[0]))
            d = self.py._data
            
            print(d.window.width)
            
            d.set_window(wbeg, wend)

            # update particle motion
            chopxy = d._chopxy()
            self.ppm.set_data(chopxy[1], chopxy[0])
            
            print(d.window.width)
            
            # update surfaces
            e = d.Py(**kwargs)
            for coll in self.sc_surf.collections: 
                self.ax2.collections.remove(coll)
            self.sc_surf = self.ax2.contourf(*e._grid, e.sc.vals, 26, cmap='magma')
            for coll in self.xc_surf.collections:
                self.ax3.collections.remove(coll)
            self.xc_surf = plt.contourf(*e._grid, e.xc.vals, 26, cmap='magma')

            
            
            
            plt.draw()
        #
        #     # if reached here process the click
        #     wbeg, wend = sorted((self.x1, self.x2))
        #
        # self.SplitWave.set_window(wbeg, wend)
   
if __name__ == "__main__":
    
    a = sw.SplitWave(split=(30,1.2), noise=0.03).Py(maxlag=4)
    explore = Explorer(a)