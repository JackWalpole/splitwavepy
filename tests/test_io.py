#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test the fundamental routines
"""

import os
import unittest, pytest
import numpy as np
import numpy.testing as npt

import splitwavepy as sw


class CoreTestCases(unittest.TestCase):
    
    def test_pickle_EigenM_io(self):
        
        # generate an eigenM object
        a = sw.EigenM()
        
        # save to file
        filename = 'temp.eigm'
        a.save(filename)
        
        # load
        b = sw.load(filename)
        
        # check they are the same
        assert a == b
        
        # cleanup
        try:
            os.remove(filename)
        except OSError:
            print("Error: cleanup failed for some reason")
        


def suite():
    asuit = unittest.makeSuite(CoreTestCases, 'test')
    return asuit
    
if __name__ == '__main__':
    unittest.main(defaultTest='suite')
