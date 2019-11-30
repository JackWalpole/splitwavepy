# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from .core.data import Data
from .core.measure import Meas
from .core.group import Group
from .core import core

# functions
from .funcs import synth, synthmany, fromst, fromsac

__all__ = ["Data, Meas, Group, core, synth, synthmany, fromst, fromsac"]



# make a commandline tool?
if __name__ == '__main__':
    # do something with commandline arguments?
    pass

