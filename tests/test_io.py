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
        a = sw.EigenM(delta=0.1)
        
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
            
    def test_pickle_TransM_io(self):
        
        # generate an eigenM object
        a = sw.TransM(delta=0.1,pol=0)
        
        # save to file
        filename = 'temp.trsm'
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

    def test_pickle_CrossM_io(self):
        
        # generate an eigenM object
        a = sw.CrossM(delta=0.1)
        
        # save to file
        filename = 'temp.crsm'
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
