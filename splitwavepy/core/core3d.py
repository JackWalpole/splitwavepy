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
    


# def pca(data):
#     """
#     Principal component analysis
#     Returns direction of strongest component in degrees anti-clockwise from x
#     """
#     w,v = np.linalg.eig(np.cov(data))
#     m = np.argmax(w)
#     return np.rad2deg(np.arctan2(v[1,m],v[0,m]))


    
