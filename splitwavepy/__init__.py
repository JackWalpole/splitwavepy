# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core.data import Data
# from .core.measure import Py
# from .core.pair import Pair
# from .core.trio import Trio
# from .measure.eigenM import EigenM
# from .measure.transM import TransM
# from .measure.crossM import CrossM
# from .measure.eig3dM import Eig3dM
# from .core import core, geom
# from .core.window import Window

__all__ = ["Data"]

# user layer functions
def synth(**kwargs):
    return Data(**kwargs).Meas(**kwargs)

def stream(**kwargs):
    # do something with ObsPy
    pass

# make a commandline tool?
if __name__ == '__main__':
    # do something with commandline arguments?
    pass

