# -*- coding: utf-8 -*-
"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ..core import core
from .eigenM import SC
from .xcorrM import XC
# from .measure import Measure

import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec


class Q(SC,XC):
    
    """
    Measure splitting using 
    
    """
    
    def __init__(self, data, **kwargs):    
        self.data = data    
        self.sc = data.SC(**kwargs)
        self.xc = data.XC(**kwargs)
        
    