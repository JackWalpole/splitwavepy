from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import .core
import .plotting
import .eigval

class Pair:
    """
    The Pair is a class to store two traces in the x and y directions.
    Methods are included to facilitate analysis on this Pair of traces.
    If data is not provided on initiation will return a ricker wavelet with noise.
    Usage: Pair()     => create Pair of synthetic data
           Pair(data) => creates Pair from two traces stored as rows in numpy array data
           Pair(x,y) => creates Pair from two traces stored in numpy arrays x and y.
    Optional:
        - set sample interval using delta=x.  Where x is a float.  Default x=1.0.
        - set orientation of trace1 using angle=x.  Where x is a float representing direction clockwise from North (or SV "up" if in ray frame).  Default is x=0.0.
    """
    def __init__(self,*args,delta=None,angle=None):
        
        if len(args) == 0:                      
            self.data = core.synth()            
        elif len(args) == 1:            
            self.data = args[0]       
        elif len(args) == 2:            
            self.data = np.vstack((args[0],args[1]))     
        else: 
            raise Exception('Unexpected number of arguments')
                    
        # some sanity checks
        if self.data.ndim != 2:
            raise Exception('data must be two dimensional')
        if self.data.shape[0] != 2:
            raise Exception('data must contain two traces in two rows')
        if self.data.shape[1]%2 == 0:
            raise Exception('traces must have odd number of samples')

        if delta is None:
            self.delta = 1.
        else:
            self.delta = float(delta)
            
        if angle is None:
            self.angle = 0.
        else:
            self.angle = float(angle)
            
    # methods
    def plot(self):
        plotting.plot_data(self.data)
    
    def split(self,degrees,tlag):
        """
        Applies splitting operator (phi,dt) to Pair.
        
        Rotates data so that trace1 is lined up with degrees (and trace2 90 degrees clockwise).
        Applies a relative time shift by the nearest even number of samples to tlag,
        trace1 is shifted tlag/2 forward in time, and trace2 tlag/2 backward in time.
        Then undoes the original rotation.
        """
        # convert time shift to nsamples -- must be even
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        # find appropriate rotation angle -- convert to counter-clockwise
        rangle = -(degrees - self.angle)
        # apply splitting
        self.data = core.split(self.data,rangle,nsamps)
        return self
    
    def unsplit(self,degrees,tlag):
        """
        Applies reverse splitting operator (phi,dt) to Pair.
        
        Rotates data so that trace1 is lined up with degrees (and trace2 90 degrees clockwise).
        Applies a relative time shift by the nearest even number of samples to tlag,
        trace1 is shifted tlag/2 backward in time, and trace2 tlag/2 forward in time.
        Then undoes the original rotation.
        """
        # convert time shift to nsamples -- must be even
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        # find appropriate rotation angle
        rangle = degrees - self.angle
        self.data = core.unsplit(self.data,rangle,nsamps)
        
    def rotateto(self,degrees):
        """
        Rotate data so that trace1 lines up with *degrees*
        """
        # find appropriate rotation angle
        rangle = degrees - self.angle
        self.data = core.rotate(self.data,rangle)
        self.angle = degrees
        
    def lag(self,tlag):
        """
        Relative shift trace1 and trace2 by tlag seconds
        """
        # convert time shift to nsamples -- must be even
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        self.data = core.lag(self.data,nsamps)
        
    def window(self,width):
        # convert width to nsamples -- must be odd
        nsamps = int(tlag / self.delta)
        nsamps = nsamps if nsamps%2==1 else nsamps + 1
        self.data = core.window(self.data,nsamps)
        
    # def grideigval(self, maxshift=None, window=None, stepang=None, stepshift=None):
    #     return eigval.grideigval(self.data)

# def time2nsamps(time,delta):
#     """
#     Convert a time to number of samples given knowledge of sample interval
#     """
#     return time/delta
#
# def nsamps2time(nsamps,delta):
#     """
#     Convert a number of samples to time given knowledge of sample interval
#     """
    