# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, core3d, geom, io
from .pair import Pair
from .window import Window

import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
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
        - window = None (default) | Window object
    
    Geometry Keyword Arguments:
        - ray = (azimuth, incidence) | tuple
        - geom = 'cart' (x,y,z) [default] | 'geo' (az,inc,r) | 'ray' (P,SH,SV) 
        - cmpvecs = np.ones(3) | custom numpy array
        - rcvloc = None
        - srcloc = None
    """
    def __init__(self,*args,**kwargs):
        
        # important to do first
        self.delta = 1.
        if ('delta' in kwargs): self.delta = kwargs['delta']
        
        # if pol specified set
        if ('pol' in kwargs): 
            self.set_pol(kwargs['pol'])
        else:
            self.set_pol()
        
        # geometry info 
        self.geom = 'geo'
        self.cmpvecs = np.eye(3)
        self.rayvecs = np.eye(3) # actual ray is assumed to be third column vector
        if ('geom' in kwargs): self.geom = kwargs['geom']
        if ('cmpvecs' in kwargs): self.cmpvecs = kwargs['cmpvecs']
        if ('ray' in kwargs): self.set_ray(*kwargs['ray'])            
                       
        # if no args make synthetic
        if len(args) == 0: 
            self.x, self.y, self.z = _synth(**kwargs)
            # rotate data to ray (but not components)

                      
        # otherwise read in data                
        elif len(args) == 3:
            if not (isinstance(args[0], np.ndarray) & 
                    isinstance(args[1], np.ndarray) &
                    isinstance(args[2], np.ndarray)):
                raise TypeError('expecting numpy arrays')         
            self.x, self.y, self.z = args[0], args[1], args[2]
        else: raise Exception('Unexpected number of arguments')
                    
        # some sanity checks
        if self.x.ndim != 1: raise Exception('data must be one dimensional')
        if self.x.size%2 == 0: raise Exception('data must have odd number of samples')
        if (self.x.size != self.y.size): raise Exception('x and y must be the same length')
                   
        # Must have a window
        self.set_window(**kwargs)
        
        # Must have a ray
        if 'ray' not in kwargs: kwargs['ray'] = (0,0)
        self.set_ray(*kwargs['ray'])
        
        # source and receiver location info
        if ('srcloc' in kwargs): self.srcloc = kwargs['srcloc']     
        if ('rcvloc' in kwargs): self.rcvloc = kwargs['rcvloc']
        if ('rayloc' in kwargs): self.rayloc = kwargs['rayloc']

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
        
    def split(self,fast,lag):
        """
        Applies splitting operator orthogonal to ray vector.

        .. warning:: shortens trace length by *lag*.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag,self.delta,mode='even')
        # find appropriate rotation angle
        origvecs = self.cmpvecs
        self.rotate2ray()
        # apply splitting
        self.x, self.y , self.z = core3d.split(self.x,self.y,self.z,fast,samps)
        self.rotateto(origvecs)

    def unsplit(self,fast,lag):
        """
        Reverses splitting operator orthogonal to ray vector.

        .. warning:: shortens trace length by *lag*.
        """
        # convert time shift to nsamples -- must be even
        samps = core.time2samps(lag,self.delta,mode='even')
        # find appropriate rotation angle
        origvecs = self.cmpvecs
        self.rotate2ray()
        # apply splitting
        self.x, self.y , self.z = core3d.unsplit(self.x,self.y,self.z,fast,samps)
        self.rotateto(origvecs)
        


    def rotate2ray(self):
        """
        Rotate data with shear plane normal to 3 axis and project 1 from "up" direction.
        """
        self.rotateto(self.rayvecs)                
    
    def rotate2eye(self):
        self.rotateto(np.eye(3))
        
    def rotate2eig(self):
        self.rotateto(self.eigvecs())
        
    def p_rotate(self):
        """
        Rotate data given P in window
        """
        p = self.eigvecs()[:,0]
        self.set_ray(p)
        self.rotate2ray()       

    def rotateto(self,vecs):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        if not np.allclose(np.eye(3),np.dot(vecs,vecs.T)):
            raise Exception('vecs must be orthogonal 3x3 matrix')
        # define the new cmpvecs
        backoff = self.cmpvecs
        self.cmpvecs = vecs
        rot = np.dot(vecs.T,backoff)
        # rotate data and ray to vecs
        xyz = np.dot(rot,self.data())
        self.x, self.y, self.z = xyz[0], xyz[1], xyz[2]
        # reset label
        self.set_labels()
        

    # def rotz(self,degs):
    #     """Rotate about z axis."""
    #     rads = math.radians(degs)
    #     rot = np.array([[ np.cos(phi), np.sin(phi), 0],
    #                     [-np.sin(phi), np.cos(phi), 0],
    #                     [           0,           0, 1]])
    #     rotateto(rot)

        
    # Set things
    
    def set_ray(self,*args):
        """
        Set the ray direction.
        
        Usage:
        
        set_ray() -- estimate ray from eigenvectors
        set_ray( np.array([ rayx, rayy, rayz]) ) -- set ray direction
        set_ray( np.array([[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]])) -- ray in column 3
        set_ray(az,inc) -- set ray from azimuth and incidence angle 
        
        if args is empty will set_ray to eigenvectors.
        """
        
        if len(args) == 0:
            self.ray = self.eigvecs()[:,2]
            sv = geom.vunit(geom.vreject([0,0,1],self.ray))
            sh = np.cross(sv,self.ray)
            self.rayvecs = np.column_stack((sv,sh,self.ray))
            return
        
        if len(args) == 1:
            if not isinstance(args[0],np.ndarray):
                raise TypeError('expecting numpy array')
            elif args[0].shape == (3,):
                self.ray = geom.vunit(args[0])
                sv = geom.vunit(geom.vreject([0,0,1],self.ray))
                sh = np.cross(sv,self.ray)
                self.rayvecs = np.column_stack((sv,sh,self.ray))
            elif args[0].shape == (3,3):
                self.rayvecs = args[0]
                self.ray = self.rayvecs[:,2]
            else:
                raise Exception('must be a shape (3,) or (3,3) array')
            return
            
        if len(args) == 2:
            az, inc = math.radians(args[0]), math.radians(args[1])
            sinc, cinc = math.sin(inc), math.cos(inc)
            saz, caz = math.sin(az), math.cos(az) 
            # Left Handed, N, E, Up
            self.rayvecs = np.array([[ cinc*caz, -saz, sinc*caz],
                                     [ cinc*saz,  caz, sinc*saz],
                                     [    -sinc,    0,     cinc]])
            self.ray = self.rayvecs[:,2]
            return
            
        # if reached here then raise exception
        raise ValueError('Unexpected arguments')
                
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
            if start > end: raise ValueError('start is larger than end')
            time_centre = (start + end)/2
            time_width = end - start
            tcs = core.time2samps(time_centre,self.delta)
            offset = tcs - self._centresamp()
            # convert time to nsamples -- must be odd
            width = core.time2samps(time_width,self.delta,'even') + 1
            self.window = Window(width,offset,**kwargs) 
            return
        else:
            raise Exception ('unexpected number of arguments')

    def set_labels(self,*args):
        if len(args) == 0:
            if np.allclose(self.cmpvecs,np.eye(3),atol=1e-02):
                if self.geom == 'geo': self.cmplabels = ['North','East','Up'] #
                # elif self.geom == 'ray': self.cmplabels = ['Vertical','Horizontal','Ray']
                # elif self.geom == 'cart': self.cmplabels = ['X','Y','Z']
                # else: self.cmplabels = ['Comp1','Comp2','Comp3']
                return
            if np.allclose(self.rayvecs,self.cmpvecs):
                self.cmplabels = ['SV','SH','P']
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


    # def geom2ray(self):
    #     """Change geometry to ray frame."""
    #
    #     if self.geom == 'ray':
    #         return
    #
    #     if self.geom == 'geo':
    #         rot = geom.geo2ray
    #         self.cmpvecs = geom.geo2ray(self.cmpvecs)
    #         self.rayvecs =
    #         self.geom = 'ray'
    #         self.set_labels()
    #
    # def geom2geo(self):
    #     """Change geometry to geographic frame."""
    #
    #     if self.geom == 'geo':
    #         return
    #
    #     if self.geom == 'ray':
    #         rot = geom.ray2geo()
    #         self.rotateto(rot)
    #         self.cmpvecs = np.eye(3)
    #         self.geom = 'geo'
    #         self.set_labels()
            
        # if self.geom ==
        #     self.
        
    
    # def rotate2cart

        
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
    noisewidth = width/4  
    
    # override defaults
    if ('pol' in kwargs): pol = kwargs['pol']   
    if ('delta' in kwargs): delta = kwargs['delta']  
    # if ('fast' in kwargs): fast = kwargs['fast']
    # if ('lag' in kwargs): lag = kwargs['lag']
    if ('split') in kwargs: split = kwargs['split']
    if ('noise' in kwargs): noise = kwargs['noise']   
    if ('nsamps' in kwargs): nsamps = kwargs['nsamps']   
    if ('width' in kwargs): width = kwargs['width'] 
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
    
    if 'ray' in kwargs:
        if not isinstance(kwargs['ray'], tuple):
            raise Exception('ray must be a tuple (azi,inc)')
        if len(kwargs['ray']) != 2:
            raise Exception('ray must be length 2 (azi,inc)')
        az, inc = math.radians(kwargs['ray'][0]), math.radians(kwargs['ray'][1])
        sinc, cinc = math.sin(inc), math.cos(inc)
        saz, caz = math.sin(az), math.cos(az)
        rot = np.array([[-cinc*caz,  saz, sinc*caz],
                        [-cinc*saz, -caz, sinc*saz],
                        [     sinc,    0,     cinc]])
        xyz = np.dot(rot,np.vstack((x,y,z)))
        x, y, z = xyz[0], xyz[1], xyz[2]
    
    return x,y,z
