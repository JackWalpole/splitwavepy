# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

class Window:
    """
    Instantiate a Window defined relative to centre of a window of flexible size.
    
    args

    - width    | nsamps length of window,
    - offset   | nsamps offset from centre of window,    
    
    kwargs
    
    - tukey   | fraction of window to cosine taper (from 0 to 1).
    """
    
    def __init__(self,width,offset=0,tukey=None):
        # ensure width is odd 
        if width%2 != 1:
            raise Exception('width must be an odd integer')
        self.width = width
        self.offset = offset
        self.tukey = tukey
    
    def start(self,samps):
        """
        Return start sample of window.
        """
        hw = int(self.width/2)
        if samps%2 != 1:
            raise Exception('samps must be odd to have definite centre')
        else:
            centre = np.int(samps/2)
            return centre + self.offset - hw

    def end(self,samps):
        """
        Return end sample of window.
        """
        hw = int(self.width/2)
        if samps%2 != 1:
            raise Exception('samps must be odd to have definite centre')
        else:
            centre = int(samps/2)
            return centre + self.offset + hw
    
    def centre(self,samps):
        """
        Return centre sample of window.
        """
        if samps%2 != 1:
            raise Exception('samps must be odd to have definite centre')
        else:
            centre = int(samps/2)
            return centre + self.offset       

    def asarray(self,samps):
                
        # sense check -- is window in range?
        if self.end(samps) > samps:
            raise Exception('Window exceeds max range')        
        if self.start(samps) < 0:
            raise Exception('Window exceeds min range')
        
        # sexy cosine taper
        if self.tukey is None:
            alpha = 0.
        else:
            alpha = self.tukey
        tukey = signal.tukey(self.width,alpha=alpha)        
        array = np.zeros(samps)
        array[self.start(samps):self.end(samps)+1] = tukey
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
        # ensure resize is even
        self.width = self.width + core.even(resize)
        
    def retukey(self,tukey):
        self.tukey = tukey
        
    # def plot(self,samps):
    #     plt.plot(self.asarray(samps))
    #     plt.show()
        
    # Comparison
    
    def __eq__(self, other) :
        if self.__class__ != other.__class__: return False
        if set(self.__dict__) != set(other.__dict__): return False
        return True
        
    # def save(self,filename):
    #     """
    #     Save just the data for future referral
    #     """
    #     io.save(self,filename)
