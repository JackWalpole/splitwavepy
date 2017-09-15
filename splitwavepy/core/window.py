from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

class Window:
    """
    Instantiate a Window 
    """
    def __init__(self,centre,width,tukey=None):
        self.centre = int(centre)
        # ensure width is odd (round up)
        width = int(width)
        self.width = width if width%2==1 else width + 1
        hw = int(self.width/2)
        self.start = self.centre - hw
        self.end = self.centre + hw
        self.tukey = tukey

    def asarray(self,nsamps=None):
        
        if nsamps is None:
            nsamps = self.start + self.end
        
        # sense check -- is window in range?
        if self.end > nsamps:
            raise Exception('Window exceeds max range')        
        if self.start < 0:
            raise Exception('Window exceeds min range')
        
        # sexy cosine taper
        if self.tukey is None:
            alpha = 0.
        else:
            alpha = self.tukey
        tukey = signal.tukey(self.width,alpha=alpha)        
        array = np.zeros(nsamps)
        array[self.start:self.end+1] = tukey
        return array
                
    def shift(self,shift):
        """
        +ve moves N samples to the right
        """
        self.centre = self.centre + int(shift)
        hw = int(self.width/2)
        self.start = self.centre - hw
        self.end = self.centre + hw
        
    def resize(self,resize):
        """
        +ve adds N samples to the window width
        """        
        # ensure resize is even (round up)
        resize = int(resize)
        resize = resize if resize%2==0 else resize + 1
        self.width = self.width + int(resize)
        hw = int(self.width/2)
        self.start = self.centre - hw
        self.end = self.centre + hw
        
    def plot(self,nsamps=None):
        plt.plot(self.asarray(nsamps))
        plt.show()