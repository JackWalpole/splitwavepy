# -*- coding: utf-8 -*-
"""
The seismic data class
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core, core3d, io
from ..core.pair import Pair
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
        
        # two or three component data
        if kwargs['ncomps'] == 2:
            base = core
        elif kwargs['ncomps'] == 3:
            base = core3d
        else:
            raise ValueError('ncomps must be 2 or 3')
            
        #     
        
            
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
                
    # Common methods
    
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
            # if not len(args[0]) == 2: raise Exception('list must be length 2')
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
        return np.vstack(( self.x, self.y, self.z))
        
    def get_pol(self):
        """Return polarisation vectors constrained normal to ray"""
        # rotate to ray
        data = self.chop()
        data.rotate2ray()
        # project to plane normal to ray
        proj = np.array([[1,0,0],
                         [0,1,0],
                         [0,0,0]])
        xyz = np.dot(proj,data.data())
        data.x, data.y, data.z = xyz[0], xyz[1], xyz[2]
        # rotate to I
        data.rotate2eye()
        # find eigvecs
        eigvecs = data.eigvecs()
        # return
        return eigvecs


    # def eigdata(self):
    #     """Return to maximum, intermediate, and minimum directions."""
    #     # rotate to I
    #     rot = self.cmpvecs.T
    #     data = self.chop().data()
    #     xyz = np.dot(rot,data)
    #     _,eigvecs = core.eigcov(xyz)
    #     rotateto(self,eigvecs)
        
    def eigvals(self):
        """Return principal component vector."""
        # rotate to I
        data = self.copy().chop()
        data.rotate2eye()
        eigvals,_ = core.eigcov(data.data())
        return(eigvals)  
        
    def eigvecs(self):
        """Return principal component vector."""
        # rotate to I
        data = self.copy().chop()
        data.rotate2eye()
        _,eigvecs = core.eigcov(data.data())
        return(eigvecs)
        

    def power(self):
        return self.x**2,self.y**2,self.z**2

    def cmpangs(self):
        """Return (az,inc) tuples in list"""
        cmp1 = self.cmpvecs[:,0]
        cmp2 = self.cmpvecs[:,1]
        cmp3 = self.cmpvecs[:,2]
        def getazi(c) : return math.degrees(math.atan2(c[1],c[0]))
        def getinc(c) : return math.degrees(math.atan2((c[0]**2+c[1]**2)**0.5,c[2]))
        return [ (getazi(cmp1),getinc(cmp1)), (getazi(cmp2),getinc(cmp2)), (getazi(cmp3),getinc(cmp3)) ]

    def chop(self):
        """
        Chop data to window
        """
        chop = self.copy()
        chop.x, chop.y, chop.z = core.chop(chop.x,chop.y,chop.z,window=chop.window)
        chop.window.offset = 0
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
        
    # plotting
    
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
        
        data = self.chop()
        data.rotate2eye()
        x, y, z = data.x, data.y, data.z
        t = data.t()
        
        # set limit
        lim = np.abs(data.data()).max() * 1.1
        if 'lims' not in kwargs: kwargs['lims'] = [-lim,lim] 
        ax.set_aspect('equal')
        ax.set_xlim(kwargs['lims'])
        ax.set_ylim(kwargs['lims'])
        ax.set_zlim(kwargs['lims'])
        
        # multi-colored
        norm = plt.Normalize(t.min(),t.max())
        points = np.array([x,y,z]).T.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = Line3DCollection(segments,cmap='plasma',norm=norm,alpha=0.7)
        lc.set_array(t)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        # plt.colorbar(line)
                
        # plot data
        # xyz = self.copy().chop().data()
        # x,y,z = xyz[0], xyz[1], xyz[2]
        # ax.plot( x, y, z)
        
        # side panel data
        ax.plot(x, y,-lim, zdir='z', alpha=0.3, color='g')
        ax.plot(x, z, lim, zdir='y', alpha=0.3, color='g')
        ax.plot(y, z, lim, zdir='x', alpha=0.3, color='g')
            
        # plot ray arrow
        rayx, rayy, rayz = self.rayvecs[0,2], self.rayvecs[1,2], self.rayvecs[2,2]
        l = 1.5*lim
        ax.quiver(0,0,0,rayx,rayy,rayz,
                  pivot='middle', color='r', length=l, alpha=0.5)
                  
        # side panel ray
        ax.quiver(lim,0,0,0,rayy,rayz,alpha=0.3,color='b',pivot='middle',length=l*math.sqrt(rayy**2+rayz**2))
        ax.quiver(0,lim,0,rayx,0,rayz,alpha=0.3,color='b',pivot='middle',length=l*math.sqrt(rayx**2+rayz**2))
        ax.quiver(0,0,-lim,rayx,rayy,0,alpha=0.3,color='b',pivot='middle',length=l*math.sqrt(rayx**2+rayy**2))
 
    
        # set labels
        if 'cmplabels' not in kwargs: kwargs['cmplabels'] = data.cmplabels
        ax.set_xlabel(kwargs['cmplabels'][0])
        ax.set_ylabel(kwargs['cmplabels'][1])
        ax.set_zlabel(kwargs['cmplabels'][2])
                
        # turn off tick annotation
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.axes.zaxis.set_ticklabels([])
        
        # flip axes
        ax.invert_xaxis()
        # ax.invert_zaxis()
        
        return
        
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