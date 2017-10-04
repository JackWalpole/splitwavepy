# -*- coding: utf-8 -*-
"""Subdivided icosahedral mesh generation"""
from __future__ import print_function
import numpy as np

# following: http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
# hierarchy:
# Icosphere -> Triangle -> Point

class IcoSphere:
    """
    Usage: IcoSphere(level)
    Maximum supported level = 8
    
    get started with:
    >>> A = IcoSphere(3)
    ... A.plot()
    """
    
    # maximum level for subdivision of the icosahedron
    maxlevel = 8
     
    def __init__(self, level):       
        if type(level) is not int:
            raise TypeError('level must be an integer')
        elif level < 0:
            raise Exception('level must be no less than 0')
        elif level > self.maxlevel:
            raise Exception('level larger than ' + str(self.maxlevel) + ' not supported')
        
        self.level = level
        self.points = []
        self.triangles = []
        self.npts = 0
        
        ################################
        # initialise level 1 icosahedron
        ################################
        
        # golden ration
        t = (1.0 + np.sqrt(5.0)) / 2.0
        
        # add vertices
        self._addPoint(np.array([-1, t, 0]))
        self._addPoint(np.array([ 1, t, 0]))
        self._addPoint(np.array([-1,-t, 0]))
        self._addPoint(np.array([ 1,-t, 0]))
        self._addPoint(np.array([ 0,-1, t]))
        self._addPoint(np.array([ 0, 1, t]))
        self._addPoint(np.array([ 0,-1,-t]))
        self._addPoint(np.array([ 0, 1,-t]))
        self._addPoint(np.array([ t, 0,-1]))
        self._addPoint(np.array([ t, 0, 1]))
        self._addPoint(np.array([-t, 0,-1]))
        self._addPoint(np.array([-t, 0, 1]))
        
        # make triangles
        tris = self.triangles
        verts = self.points
        # 5 faces around point 0
        tris.append(Triangle([ verts[0],verts[11], verts[5]]))
        tris.append(Triangle([ verts[0], verts[5], verts[1]]))
        tris.append(Triangle([ verts[0], verts[1], verts[7]]))
        tris.append(Triangle([ verts[0], verts[7],verts[10]]))
        tris.append(Triangle([ verts[0],verts[10],verts[11]]))
        # 5 adjacent faces
        tris.append(Triangle([ verts[1], verts[5], verts[9]]))
        tris.append(Triangle([ verts[5],verts[11], verts[4]]))
        tris.append(Triangle([verts[11],verts[10], verts[2]]))
        tris.append(Triangle([verts[10], verts[7], verts[6]]))
        tris.append(Triangle([ verts[7], verts[1], verts[8]]))
        # 5 faces around point 3
        tris.append(Triangle([ verts[3], verts[9], verts[4]]))
        tris.append(Triangle([ verts[3], verts[4], verts[2]]))
        tris.append(Triangle([ verts[3], verts[2], verts[6]]))
        tris.append(Triangle([ verts[3], verts[6], verts[8]]))
        tris.append(Triangle([ verts[3], verts[8], verts[9]]))
        # 5 adjacent faces
        tris.append(Triangle([ verts[4], verts[9], verts[5]]))
        tris.append(Triangle([ verts[2], verts[4],verts[11]]))
        tris.append(Triangle([ verts[6], verts[2],verts[10]]))
        tris.append(Triangle([ verts[8], verts[6], verts[7]]))
        tris.append(Triangle([ verts[9], verts[8], verts[1]]))
        
        ########################################
        # refine triangles to desired mesh level
        ########################################
        
        for l in range(self.level):
            midPointDict = {}
            faces = []
            for tri in self.triangles:
                # replace triangle by 4 triangles
                p = tri.pts
                a = self._getMiddlePoint(p[0], p[1], midPointDict)
                b = self._getMiddlePoint(p[1], p[2], midPointDict)
                c = self._getMiddlePoint(p[2], p[0], midPointDict)
                faces.append(Triangle([p[0], a, c]))
                faces.append(Triangle([p[1], b, a]))
                faces.append(Triangle([p[2], c, b]))
                faces.append(Triangle([a, b, c]))
            # once looped thru all triangles overwrite self.triangles
            self.triangles = faces
        
        self.nfaces = len(self.triangles)
            
        # check that npts and nfaces are as expected
        expected_npts = calculate_npts(self.level)
        expected_nfaces = calculate_nfaces(self.level)
        if self.npts != calculate_npts(self.level):
            raise Exception('npts '+str(self.npts)+' not as expected '+str(expected_npts))
        elif self.nfaces != calculate_nfaces(self.level):
            raise Exception('nfaces '+str(self.nfaces)+' not as expected '+str(expected_nfaces))
    
    def _addPoint(self, xyz):
        """Add point to self.points"""
        self.points.append(Point(self.npts, xyz))
        self.npts += 1
        
    def _getMiddlePoint(self, p1, p2, midPointDict):
        """return Point"""
        if not isinstance(p1, Point) or not isinstance(p2, Point):
            raise TypeError('p1 and p2 must be Points')        
        # does point already exist?
        key = tuple(sorted([p1.idx, p2.idx]))
        if key in midPointDict:
            # point exists
            pass
        else:
            # point is new
            self._addPoint((p1.xyz + p2.xyz)/2)
            midPointDict[key] = self.points[-1] 
        return midPointDict[key]
        
    def plot(self):
        """Matplotlib 3D plot of mesh"""
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xyz = np.asarray([ pt.xyz for pt in self.points ])
        x = xyz[:,0]
        y = xyz[:,1]
        z = xyz[:,2]
        ts = np.asarray([ [ p.idx for p in t.pts ] for t in self.triangles ])
        ax.plot_trisurf(x,y,ts,z)
        plt.show()
        
        
        
class Triangle:
    """A triangle adjoining three adjacent points"""
    def __init__(self, pts):              
        if not isinstance(pts, list):
            raise TypeError('pts must be a list')
        elif len(pts) !=3:
            raise Exception('pts must be of length 3')    
        else:
            self.pts = pts

class Point:
    """A 3D point on the mesh"""
    def __init__(self, idx, xyz):         
        
        if type(idx) is not int:
            raise TypeError('idx must be an integer')              
        elif not isinstance(xyz,np.ndarray):
            raise TypeError('xyz must be a numpy array')
        elif xyz.size != 3:
            raise Exception('xyz must be of size 3')
        else:
            # ensure length equals 1 and add to list of points
            self.xyz = (xyz/np.linalg.norm(xyz))
            self.idx = idx
            
        
def calculate_npts(level):
    n = 2**level
    return 2 + 10 * n**2
    
def calculate_nfaces(level):
    n = 2**level
    return 20 * n**2
    
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
    x = r * np.cos(lon) * np.cos(lat)
    y = r * np.sin(lon) * np.cos(lat)
    z = r * np.sin(lat)
    return x, y, z

# def xyzToLatLonR(xyz):
#     trans = np.array([np.])
      
