# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core.data import Data
from .core.measure import Meas
from .core.group import Group

from obspy import read

# user layer functions
def synth(**kwargs):
    """Generate a synthetic measurement."""
    return Data(**kwargs).Meas(**kwargs)
    
def synthmany(n, **kwargs):
    """Generate n synthetic measurements."""
    return sw.Group([ sw.synth(**kwargs) for x in range(n) ])
    

    

def stream(**kwargs):
    # do something with ObsPy
    pass