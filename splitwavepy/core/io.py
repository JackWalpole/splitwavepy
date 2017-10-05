# -*- coding: utf-8 -*-
"""
Save and load things
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pickle

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
    Load an EigenM object
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)