# -*- coding: utf-8 -*-
"""
Ray geometry stuff
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# import taup from obsp for ray calculation
import numpy as np
import math


# classes

class Ray:
    """
    collection of attributes for a raypath
    """
    # self.srcloc =
    # self.rcvloc =
    # self.path =

class Point:
    """
    numpy xyz vector
    """
    
    def __init__(self,xyz):
        self.xyz = xyz

# useful co-ordinate transforms

def cart2geo(x, y, z):
    """convert x y z cartesian coordinates to latitude longitude radius
       xyz is a numpy array, a right handed co-ordinate system is assumed with 
       -- x-axis going through the equator at 0 degrees longitude
       -- y-axis going through the equator at 90 degrees longitude 
       -- z-axis going through the north pole."""
    r = np.sqrt(x**2 + y**2 + z**2)
    lon = np.rad2deg(np.arctan2(y,x))
    lat = np.rad2deg(np.arcsin(z/r))
    return lat, lon, r

def geo2cart(lat, lon, r):
    """convert latitude longitude radius to x y z cartesian coordinates
       xyz is a numpy array, a right handed co-ordinate system is assumed with 
       -- x-axis going through the equator at 0 degrees longitude
       -- y-axis going through the equator at 90 degrees longitude 
       -- z-axis going through the north pole."""
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    x = r * np.cos(lon) * np.cos(lat)
    y = r * np.sin(lon) * np.cos(lat)
    z = r * np.sin(lat)
    return x, y, z


    
# def sph2cart(r,az,inc):
#     """
#     convert from r,az,inc to local x,y,z
#     azimuth measured clockwise from local North, incidence angle measured from vertical down direction
#     """
#     az = np.deg2rad(az)
#     inc = np.deg2rad(inc)
#     x = r * np.cos(az) * np.sin(inc)
#     y = r * np.sin(az) * np.sin(inc)
#     z = r * np.cos(inc)
#     return x, y, z
#
#
# def cart2sph(x,y,z):
#     XsqPlusYsq = x**2 + y**2
#     r = np.sqrt(XsqPlusYsq + z**2)
#     az = np.arctan2(y,x)
#     inc = np.arctan2(np.sqrt(XsqPlusYsq),z)
#     az = np.rad2deg(az)
#     inc = np.rad2deg(inc)
#     return r, az, inc
    
# def ray2xyz(v,reverse=False):
#     """
#     conversion of vector v from ray co-ordinates to x,y,z
#     ray co-ordinates form a right-handed orthogonal system in the P, SV, and SH (1,2,3) directions
#     xyz co-ordinates form a right-handed orthogonal system in the principal geographic directions
#     for reverse transform use reverse=True
#     """
#
# def ray2any(v,reverse=False):
#     """
#     conversion of vector v from ray co-ordinates to user defined co-ordinate system
#     ray co-ordinates form a right-handed orthogonal system in the P, SV, and SH (1,2,3) directions
#     user defined co-ordinates must be provided in a 3x3 matrix with columns corresponding to directions of interest (e.g. slab frame: dip, strike, and slab normal)
#     for reverse transform use reverse=True
#     """
#
# def xyz2any(v,B,reverse=False):
#     """
#     conversion of vector v from xyz to user defined co-ordinate system
#     xyz co-ordinates form a right-handed orthogonal system in the principal geographic directions
#     user defined co-ordinates must be provided in a 3x3 matrix with columns corresponding to directions of interest (e.g. slab frame: dip, strike, and slab normal)
#     for reverse transform use reverse=True
#     """
#


#

def enu2psv(az,inc):
    """return matrix to convert local ENU coordinates to PSV coordinates
       given local azimuth and inclination of ray"""
    az = np.deg2rad(az)
    saz = np.sin(az)
    caz = np.cos(az)
    inc = np.deg2rad(inc)
    sinc = np.sin(inc)
    cinc = np.cos(inc)

    m = np.array([[sinc*saz, -caz, -cinc*saz],
                  [sinc*caz,  saz, -cinc*caz],
                  [    cinc,    0,      sinc]])
    return m.T

def psv2enu(az,inc):
    """return matrix to convert PSV coordinates to local ENU coordinates 
       given local azimuth and inclination of ray"""
    az = math.radians(az)
    inc = math.radians(inc)
    sinc = math.sin(inc)
    cinc = math.cos(inc)
    saz = math.sin(az)
    caz = math.cos(az)
    m = np.array([[sinc*saz, -caz, -cinc*saz],
                  [sinc*caz,  saz, -cinc*caz],
                  [    cinc,    0,      sinc]])
    return m

# useful directions in xyz

def vnpole():
    """
    North pole
    """
    return np.array([0,0,1])
    
def vup(lat,lon):
    """
    return unit vector pointing up
    """
    return np.asarray(geo2cart(lat,lon,1))

def vnorth(lat,lon):
    """
    return unit vector pointing north
    """
    npole = vnpole()
    up = vup(lat,lon)
    return vunit(vrejection(npole,up))

def vray(lat,lon,azi,inc):
    """
    return unit vector pointing in ray direction
    """
    up = vup(lat,lon)
    north = vnorth(lat,lon)
    azi = -np.deg2rad(azi)
    # convert north to azi vector
    azi = np.dot(rotation_matrix(up,azi),north)
    # horizontal transverse
    hor = np.cross(up,azi)
    # THE RAY
    # convert incidence angle to vector
    inc = np.deg2rad(inc)
    inc = np.dot(rotation_matrix(hor,inc),up)
    return inc
    
    
# def ray_outgoing(lat,lon,depth,az,tkoff):


#
# def ray(point1,point2):
#     """
#     return unit vector pointing from point 1 to point2
#     """
#
#
# def trans_theta(point,ray,theta):
#     """
#     return vector at point p, transverse to ray, and at angle theta degrees clockwise from  "up"
#     useful angles include: up (0), horizontal (+/-90), polarisation, fast, slow
#     if ray is "up"
#     """


# common rotations

def rz(a,phi):
    # rotate vector a phi radians clockwise around the z axis (if z pointing down, x->N, y->E)
    rot=np.array([[np.cos(phi), -np.sin(phi), 0],
                [np.sin(phi), np.cos(phi), 0],
                [0,0,1]])
    return np.dot(rot,a)

def ry(a,phi):
    # rotate vector a phi radians clockwise around the y axis (if z pointing down, x->N, y->E)
    rot=np.array([[np.cos(phi), 0, np.sin(phi)],
                [0,1,0],
                [-np.sin(phi), 0, np.cos(phi)]])
    return np.dot(rot,a)

def rx(a,phi):
    # rotate vector a phi radians clockwise around the y axis
    rot=np.array([[1,0,0],
                [0, np.cos(phi), -np.sin(phi)],               
                [0, np.sin(phi), np.cos(phi)]])
    return np.dot(rot,a)
    
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    x rotated pi/2 radians about the z axis points in the +ve y-axis direction (right handed system)
    downloaded from: https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

# vector basics
# works with numpy arrays

# def vangle(a,b):
#     # Return the angle between two vectors
#     return np.arccos(round(np.dot(a,b) / (np.linalg.norm(a)*np.linalg.norm(b)),15))

def vangle(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793

        Downloaded from:https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
    """
    v1_u = vunit(v1)
    v2_u = vunit(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def vunit(a):
    # Return a unit vector in the direction of a
    return a/np.linalg.norm(a)

def vreject(a,b):
    # The vector rejection of a on b (the bit of a perpendicular to b)
    # Equivalent to projection onto the plane normal to b
    return a - (np.linalg.norm(a)*np.cos(vangle(a,b))) * vunit(b)
    
# fast direction stuff
 
def phiray2geo(phi,az,inc):
    """
    convert phi in ray frame to phi in geographic frame
    phi is angle measured in ray frame
    az, inc is azimuth and inclination of ray at station
    """
    
    # get fast axis in ray co-ordinates
    phi = np.deg2rad(phi)
    fast = [0, -np.sin(phi), np.cos(phi)]
    # convert fast to geographic co-ordinates
    ray2geo = psv2enu(az,inc)
    fast = np.dot(ray2geo,fast)    
    # measure angle from east and north parts
    ang = np.rad2deg(np.arctan2(fast[0],fast[1]))
    return (ang+3690)%180-90



def phigeo2ray(phi,az,inc):
    """
    convert phi in ray frame to phi in geographic frame
    phi is angle measured in ray frame
    az, inc is azimuth and inclination of ray at station
    """
    
    # find where up vector at point related to phi in geographic co-ordinates intersects the plane
    phi = np.deg2rad(phi)
    fastfloor = [np.sin(phi),np.cos(phi),0]
    ray2geo = psv2enu(az,inc)
    geo2ray = ray2geo.T
    raygeo = np.dot(ray2geo,[1,0,0])
    up = -(raygeo[0]*fastfloor[0] + raygeo[1]*fastfloor[1])/raygeo[2]
    fastgeo = [fastfloor[0],fastfloor[1],up]
    fastray = np.dot(geo2ray,fastgeo)
    ang = np.rad2deg(np.arctan2(-fastray[1],fastray[2]))
    return (ang+3690)%180-90
