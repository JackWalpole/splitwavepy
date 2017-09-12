import core
import plotting
import eigval

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
    """
    def __init__(self,*args,delta=None):
        
        if len(args) == 0:                      
            self.data = synth()            
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
            
    # methods
    def plot(self):
        plotting.plot_data(self.data)
    
    def split(self,degrees,nsamps):
        self.data = core.split(self.data,degrees,nsamps)
        return self
    
    def unsplit(self,degrees,nsamps):
        self.data = core.unsplit(self.data,degrees,nsamps)
        return self
        
    def rotate(self,degrees):
        self.data = core.rotate(self.data,degrees)
        return self
        
    def lag(self,tlag):
        nsamps = int(tlag / self.delta)
        # must be even
        nsamps = nsamps if nsamps%2==0 else nsamps + 1
        self.data = core.lag(self.data,nsamps)
        return self
        
    def window(self,width):
        nsamps = int(tlag /self.data)
        # must be odd
        nsamps = nsamps if nsamps%2==1 else nsamps + 1
        self.data = core.window(self.data,nsamps)
        return self
        
    def grideigval(self, maxshift=None, window=None, stepang=None, stepshift=None):
        return eigval.grideigval(self.data)