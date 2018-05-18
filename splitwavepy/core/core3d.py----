# -*- coding: utf-8 -*-
"""
Routines to work with 3-component data
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .window import Window
from . import geom
from . import core

import numpy as np
from scipy import signal
import math

def lag(x,y,z,nsamps):
    """
    Lag t1 nsamps/2 to the left and
    lag t2 nsamps/2 to the right.
    nsamps must be even.
    If nsamps is negative t1 shifted to the right and t2 to the left.
    This process truncates trace length by nsamps and preserves centrality.
    Therefore windowing must be used after this process 
    to ensure even trace lengths when measuring splitting.
    """
    if nsamps == 0:
        return x,y,z
    elif nsamps%2 != 0:
        raise Exception('nsamps must be even')
    
    x, y = core.lag(x,y,nsamps)
    z = z[int(abs(nsamps)/2):-int(abs(nsamps)/2)]
    
    return x, y, z
        
def rotate(x,y,z,degrees):
    """row 0 is x-axis and row 1 is y-axis,
       rotates from x to y axis
       e.g. N to E if row 0 is N cmp and row1 is E cmp"""
    ang = math.radians(degrees)
    rot = np.array([[ math.cos(ang), math.sin(ang), 0],
                    [-math.sin(ang), math.cos(ang), 0],
                    [             0,             0, 1]])
    xyz = np.dot(rot,np.vstack((x,y,z)))
    return xyz[0], xyz[1], xyz[2]

def split(x,y,z,degrees,nsamps):
    """Apply forward splitting and rotate back"""
    x,y,z = rotate(x,y,z,degrees)
    x,y,z = lag(x,y,z,nsamps)
    x,y,z = rotate(x,y,z,-degrees)
    return x,y,z

def unsplit(x,y,z,degrees,nsamps):
    """Apply inverse splitting and rotate back"""
    return split(x,y,z,degrees,-nsamps)
    
def chop(*args,**kwargs):
    """Chop trace, or traces, using window"""
    return core.chop(*args,**kwargs)
    
## Measurement 
   
def eigcov(data):
    """
    Return eigen values and vectors of covariance matrix
    """
    eigenValues, eigenVectors = np.linalg.eig(np.cov(data))
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return eigenValues, eigenVectors
        
def eigvalcov(x,y,z):
    """
    return sorted eigenvalues of covariance matrix
    lambda2 first, lambda1 second
    """
    data = np.vstack((x,y,z))
    return np.sort(np.linalg.eigvalsh(np.cov(data)))
  
def transenergy(x,y,z):
    """
    return energy
    lambda1 first, lambda2 second
    """
    energy = lambda x: np.sum(x**2)
    return energy(x), energy(y), energy(z) 
    
def synth(**kwargs):
    """return ricker wavelet synthetic data"""
    
    # defaults    
    if 'pol' not in kwargs: kwargs['pol'] = 0.
    if 'delta' not in kwargs: kwargs['delta'] = 1.
    if 'split' not in kwargs: kwargs['split'] = []
    if 'noise' not in kwargs: kwargs['noise'] = 0.001
    if 'nsamps' not in kwargs: kwargs['nsamps'] = 1001
    if 'width' not in kwargs: kwargs['width'] = 32
    if 'noisewidth' not in kwargs: kwargs['noisewidth'] = kwargs['width']/4

    # initiate data with clean ricker wavelet
    nsamps = int(kwargs['nsamps'])  
    x = signal.ricker(nsamps, kwargs['width'])
    y = np.zeros(nsamps)
    
    # rotate to polarisation 
    # negative because we are doing the active rotation of data, whereas
    # core.rotate does the passive transormation of the co-ordinate system
    x,y = core.rotate(x,y,-kwargs['pol'])

    if isinstance(kwargs['split'],tuple):
        fast, lag = kwargs['split']
        # add any splitting -- lag samples must be even
        slag = core.time2samps(lag,kwargs['delta'],mode='even')
        x,y = core.split(x,y,fast,slag)
    elif isinstance(split,list):        
        for parms in kwargs['split']:
            fast, lag = parms
            # add any splitting -- lag samples must be even
            slag = core.time2samps(lag,kwargs['delta'],mode='even')
            x,y = core.split(x,y,fast,slag)
    
    # add noise - do this last to avoid splitting the noise
    x = x + core.noise(x.size,kwargs['noise'],int(kwargs['noisewidth']))    
    y = y + core.noise(x.size,kwargs['noise'],int(kwargs['noisewidth']))
    z = core.noise(x.size,kwargs['noise'],int(kwargs['noisewidth']))
    
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
    
    return x, y, z

    
