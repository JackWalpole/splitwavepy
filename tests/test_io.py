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
        
        #make synthetic
        synth = sw.Data(delta=0.1,split=(30,1.4),noise=0.03)
        
        # generate an eigenM object
        a = synth.EigenM()
        
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
        
        #make synthetic
        synth = sw.Data(delta=0.1,split=(30,1.4),noise=0.03)
        
        # generate an TransM object
        a = synth.TransM(pol=0)
        
        # save to file
        filename = 'temp.trnm'
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

    def test_pickle_XcorrM_io(self):
        
        #make synthetic
        synth = sw.Data(delta=0.1,split=(30,1.4),noise=0.03)
        
        # generate an XcorrM object
        a = synth.XcorrM()
        
        # save to file
        filename = 'temp.xcrm'
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
