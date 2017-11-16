# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, core3d, geom, io
# from .pair import Pair
from .data import Data, WindowPicker
from .window import Window

import numpy as np
import math
from scipy import signal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import gridspec

class Trio(Data):
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
        
        # Derive from Data
        Data.__init__(self, *args, **kwargs)   
                       
        # if no args make synthetic
        if len(args) == 0: 
            self.x, self.y, self.z = core3d.synth(**kwargs)                      
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
        if (self.x.size != self.y.size or self.x.size != self.z.size): 
            raise Exception('x, y, and z must be the same length')
        
        # geometry info 
        self.geom = 'geo'
        self.cmpvecs = np.eye(3)
        self.rayvecs = np.eye(3) # actual ray is assumed to be third column vector
        if ('geom' in kwargs): self.geom = kwargs['geom']
        if ('cmpvecs' in kwargs): self.cmpvecs = kwargs['cmpvecs']
        if ('ray' in kwargs): self.set_ray(*kwargs['ray'])
                   
        # Must have a window
        self.set_window(**kwargs)
        
        # Must have a ray
        if 'ray' not in kwargs: kwargs['ray'] = (0,0)
        self.set_ray(*kwargs['ray'])
        
        # Must have a pol
        if ('pol' in kwargs): 
            self.set_pol(kwargs['pol'])
        else:
            self.set_pol()
        
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
            # if not len(args[0]) == 3: raise Exception('list must be length 3')
            if not (isinstance(args[0][0],str) and 
                    isinstance(args[0][1],str) and
                    isinstance(args[0][2],str)):
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
