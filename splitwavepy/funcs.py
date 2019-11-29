# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core.data import Data
from .core.measure import Meas
from .core.group import Group

# user layer functions
def synth(**kwargs):
    return Data(**kwargs).Meas(**kwargs)

def stream(**kwargs):
    # do something with ObsPy
    pass