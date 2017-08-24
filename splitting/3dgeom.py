"""
Ray geometry stuff
"""

# import taup from obsp for ray calculation



# classes

class Ray:
    """
    collection of attributes for a raypath
    """
    self.srcloc =
    self.rcvloc =
    self.path = 

class Point:
    """
    numpy xyz vector
    """
    
    def __init__(self,xyz):
        self.xyz = xyz

# useful co-ordinate transforms

def cart2geo(x,y,z,reverse=False):
    """
    convert from x,y,z to lat,lon,r
    convention:
    x =>   0 E,  0 N 
    y =>  90 E,  0 N
    z => NaN E, 90 N
    for reverse transform use reverse=True and give lat,lon,r in place of x,y,z
    """

def xyz2azinc(x,y,z,reverse=False):
    """
    convert from global x,y,z to azimuth, and incidence angle at Point
    azimuth mesured clockwise from local North, incidence angle measured from vertical down direction
    for reverse transform use reverse=True and give az,inc,Point in place of x,y,z
    """

def ray2xyz(v,reverse=False):
    """
    conversion of vector v from ray co-ordinates to x,y,z
    ray co-ordinates form a right-handed orthogonal system in the P, SV, and SH (1,2,3) directions
    xyz co-ordinates form a right-handed orthogonal system in the principal geographic directions
    for reverse transform use reverse=True
    """

def ray2any(v,reverse=False):
    """
    conversion of vector v from ray co-ordinates to user defined co-ordinate system
    ray co-ordinates form a right-handed orthogonal system in the P, SV, and SH (1,2,3) directions
    user defined co-ordinates must be provided in a 3x3 matrix with columns corresponding to directions of interest (e.g. slab frame: dip, strike, and slab normal)
    for reverse transform use reverse=True
    """    

def xyz2any(v,B,reverse=False):
    """
    conversion of vector v from xyz to user defined co-ordinate system
    xyz co-ordinates form a right-handed orthogonal system in the principal geographic directions
    user defined co-ordinates must be provided in a 3x3 matrix with columns corresponding to directions of interest (e.g. slab frame: dip, strike, and slab normal)
    for reverse transform use reverse=True
    """



# useful directions in xyz
    
def north(point):
    """
    return unit vector pointing north
    """
    
def up():
    """
    return unit vector pointing up
    """
   
    
def ray(point1,point2):
    """
    return unit vector pointing from point 1 to point2
    """
    

def trans_theta(point,ray,theta):
    """
    return vector at point p, transverse to ray, and at angle theta degrees clockwise from  "up"
    useful angles include: up (0), horizontal (+/-90), polarisation, fast, slow
    if ray is "up"
    """


    


# vector basics
def vunit():
    """
    return unit vector in direction v
    """
    
def vreject(a,b):
    

def 