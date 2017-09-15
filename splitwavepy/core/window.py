from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

class Window:
    """
    Instantiate a Window defined relative to centre of a window of flexible size.
    """
    def __init__(self,width,offset=0,tukey=None):
        # ensure width is odd (round up)
        width = int(width)
        self.width = width if width%2==1 else width + 1
        self.offset = offset
        self.tukey = tukey
    
    def start(self,nsamps):
        """
        Return start sample of window.
        """
        hw = int(self.width/2)
        if nsamps%2 != 1:
            raise Exception('nsamps must be odd to have definite centre')
        else:
            centre = int(nsamps/2)
            return centre + self.offset - hw

    def end(self,nsamps):
        """
        Return end sample of window.
        """
        hw = int(self.width/2)
        if nsamps%2 != 1:
            raise Exception('nsamps must be odd to have definite centre')
        else:
            centre = int(nsamps/2)
            return centre + self.offset + hw
    
    def centre(self,nsamps):
        """
        Return centre sample of window.
        """
        if nsamps%2 != 1:
            raise Exception('nsamps must be odd to have definite centre')
        else:
            centre = int(nsamps/2)
            return centre + self.offset       

    def asarray(self,nsamps):
                
        # sense check -- is window in range?
        if self.end(nsamps) > nsamps:
            raise Exception('Window exceeds max range')        
        if self.start(nsamps) < 0:
            raise Exception('Window exceeds min range')
        
        # sexy cosine taper
        if self.tukey is None:
            alpha = 0.
        else:
            alpha = self.tukey
        tukey = signal.tukey(self.width,alpha=alpha)        
        array = np.zeros(nsamps)
        array[self.start(nsamps):self.end(nsamps)+1] = tukey
        return array
                
    def shift(self,shift):
        """
        +ve moves N samples to the right
        """
        self.offset = self.offset + int(shift)
        
    def resize(self,resize):
        """
        +ve adds N samples to the window width
        """        
        # ensure resize is even (round up)
        resize = int(resize)
        resize = resize if resize%2==0 else resize + 1
        self.width = self.width + int(resize)
        
    def retukey(self,tukey):
        self.tukey = tukey
        
    def plot(self,nsamps):
        plt.plot(self.asarray(nsamps))
        plt.show()