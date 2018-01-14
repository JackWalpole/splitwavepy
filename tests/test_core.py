#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test the fundamental routines
"""

import unittest, pytest
import numpy as np
import numpy.testing as npt

import splitwavepy as sw


class CoreTestCases(unittest.TestCase):
    
    def test_lag(self):
        """lag should shift arrays relative to one another"""
        x = np.array([0,1,2])
        y = np.array([3,4,5])
        npt.assert_array_equal(sw.core.lag(x,y,2), np.array([[2],[3]]))
        npt.assert_array_equal(sw.core.lag(x,y,-2), np.array([[0],[5]]))
        
    def test_rotate(self):
        """rotation of traces (co-ordinate frame)"""
        x = np.array([0,1,1])
        y = np.array([1,0,1])
        npt.assert_array_almost_equal(sw.core.rotate(x,y,90),np.array([[1,0,1],[0,-1,-1]]))
        
    def test_round_to_int(self):
        """get nearest integer"""
        npt.assert_array_equal(np.array([0, 0, -6, 12]), 
            sw.core.near(np.array([-.1, .34, -5.78, 11.88])))
        npt.assert_array_equal(np.array([0, 0, 10, -12]),
            sw.core.even(np.array([-0.9, 0.5, 9.2, -11.1])))
        npt.assert_array_equal(np.array([1, -1, 9, -11]),
            sw.core.odd(np.array([1.1, -1.1, 8.1, -10.3])))
                
    def test_time2samps(self):
        """convert time to samples"""
        assert sw.core.time2samps(1.3, 0.1) == 13
        assert sw.core.samps2time(13, 0.1) == 1.3
                
    # def test_eigcov(self):
    #
    # def test_eigvalcov(self):
    #
    # def test_transenergy(self):
    #
    # def test_crosscorr(self):
    #
    # def test_crossconv(self):
    #
    # def test_misfit(self):
    #
    # def test_crossconvmf(self):
    #
    # def test_splittingintensity(self):
    #
    # def test_ndf(self):
    #
    # def test_ftest(self):
    #
    # def test_Q(self):
    #
    # def test_snr(self):
    #
    # def test_synth(self):
    #
    # def test_noise(self):
    #
    # def test_resample_noise(self):
    #
    # def test_min_idx(self):
    #
    # def test_max_idx(self):
        
    
def suite():
    asuit = unittest.makeSuite(CoreTestCases, 'test')
    return asuit
    
if __name__ == '__main__':
    unittest.main(defaultTest='suite')
