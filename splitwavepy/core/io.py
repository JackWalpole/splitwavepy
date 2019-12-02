# -*- coding: utf-8 -*-
"""
Save and load things
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pickle
import copy

# Saving

def save(self,filename):
    """
    Save me to a file
    """       
    with open(filename, 'wb') as f:
        pickle.dump(self,f)

# Loading

def load(filename):
    """
    Load an SC object
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)
        
# Copying

copy = copy.deepcopy
#
# # Decorator
# def respawn(func):
#     print('Here')
#     def wrapper(self, **kwargs):
#         return func(copy(self), **kwargs)
#     return wrapper

# Would be good to incorporte this info in files:
# version tracking ... e.g..
# version = None
# git = None
