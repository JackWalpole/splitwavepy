# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core.data import Data
from .core.measure import Meas
from .core.group import Group

import obspy

# user layer functions
def synth(**kwargs):
    """Generate a synthetic measurement."""
    return Data(**kwargs).Meas(**kwargs)
    
def synthmany(n, **kwargs):
    """Generate n synthetic measurements."""
    return Group([ sw.synth(**kwargs) for x in range(n) ])
    
def fromst(st, **kwargs):
    """Generate from obspy stream object."""    
    delta = st[0].stats.delta
    x, y = st[1].data, st[0].data
    return Data(x, y, delta=delta, **kwargs).Meas(**kwargs)

def fromsac(*args, **kwargs):
    st = obspy.read(*args, **kwargs)
    return fromst(st)