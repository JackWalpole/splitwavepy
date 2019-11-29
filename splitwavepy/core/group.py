# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core, io
from .data import Data
from .measure import Meas

import numpy as np
import matplotlib.pyplot as plt

# useful values
sig1 = 1 - 0.6827
sig2 = 1 - 0.9545
sig3 = 1 - 0.9973

class Group:

    def __init__(self, listM, **kwargs):
        """
        A collection of measurements (or data).
        
        args
        -----
        
        listM      list, a list of Q objects (all on same grid)
        
        keyword arguments
        -----------------
        
        weights     numpy array, user-defined weights (order and size must match list)
        
               
        example
        --------
      
        >>> S = Group(listM)           # initiate a Stack of EigenM objects
        >>> l1l2stack = S.stack()           # Lam1 / Lam2 stack
        >>> WSstack = S.wolfe_silver()      # Stack using Wolfe and Silver method
        >>> RHstack = S.restivo_helffrich() # Stack using Restivo and Helffrich method
        """
            
        # get info
        self.listM = listM
        m = self.listM[0]
        self._degs = m._degs
        self._lags = m._lags
        self._grid = m._grid
        self.units = m.units
        # get result
        l = np.exp(self.mean())
        l = l / l.max() # normalise assuming we know the answer!
        self.likelihood = l
        self.fast, self.lag = m._fast_lag_maxloc(l)
        self.dfast, self.dlag = m._errorbars(l, alpha=sig1)

        
        # checks at the border guard
        # self._deltacheck()
        # self._degscheck()
        
        # self._delta
        
        # try:
        #     self._lagscheck()
        # except:
        #     self._trimlags()
        
        
    def wolfe_silver(self, **kwargs):
        """
        Return stack using method of Wolfe and Silver (1998).

        This method stacks the lambda2 surfaces,
        pre-normalises surfaces so that minimum lambda2 = 1.
        """

        listS = [ m.lam2 / np.min(m.lam2) for m in self.listM ]
        return np.mean(listS, axis=0)

    def restivo_helffrich(self, **kwargs):
        """
        Return stack using method of Restivo and Helffrich (1999).

        This method is similar to the Wolfe and Silver method except
        error surfaces are weighted by their signal to noise ratio.
        """

        listS = [ M.lam2 / np.min(M.lam2) for M in self.listM ]

        # weight by signal to noise ratio
        weights = np.asarray([ M.data.snrRH() for M in self.listM ])

        # should apply sigmoid (?) function to weights with min to max ranging from 1 to 21 to be consistent with original paper.  Note: sheba does not bother with this.

        # if 'baz' in kwargs and kwargs['baz'] is True:
        #     # weight by backazimuthal density coverage
        #     raise Exception('not yet supported')

        # if user specified weights then apply these on top:
        if 'weights' in kwargs:
            weights = weights * kwargs['weights']

        return np.mean(listS, weights=weights)
        
    @property
    def stack(self):
        stack = np.stack([ np.log(a.likelihood) for a in self.listM])
        return stack
        
    # slower functions
    
    def sum(self):
        return np.sum(self.stack, axis=0)
       
    def mean(self):
        return np.mean(self.stack, axis=0)
        
    def std(self):
        return np.std(self.stack, axis=0)
        
    def pddf(self):
        """Pandas DataFrame Textual Results."""
        pass
        
    def animate(self):
        """video of surfaces."""
        pass
        
    def plot(self):
        pass
            
    def cluster(self):
        pass
    
    def _checktype(self):
        return all(x==type(Meas) for x in self.listM)
    

    def _gridcheck(self):
        if not core.core.check_list_same([ m._grid for m in self.listM ]):
            raise Exception('Grids not all the same.')
        else:
            return True
            
    def _deltacheck(self):
        if not core.check_list_same([ m.data._delta for m in self.listM ]):
            raise Exception('Grids not all the same.')
            
    def _lagscheck(self):
        if not core.check_list_same([ m._lags for m in self.listM ]):
            raise Exception('Grids not all the same.')
            
    def _degscheck(self):
        if not core.check_list_same([ m._degs for m in self.listM ]):
            raise Exception('Grids not all the same.')
            
    def plot(self, **kwargs):            
        fig, ax = plt.subplots()
        vals = np.flipud(self.likelihood.T)
        (b, t), (l, r) = self._degs[[0,-1]], self._lags[[0, -1]]
        ax.imshow(vals, extent=(l, r, b, t), aspect='auto', **kwargs)
        ax.set_xlim([l, r])
        ax.set_ylim([b, t])
        ax.contour(*self._grid, self.likelihood, levels=[0.05], colors='white')
        ax.set_ylabel(r'Fast Direction ($^\circ$)')
        ax.set_xlabel('Delay Time (' + self.units + ')')
        # marker
        fast, dfast, lag, dlag = self.fast, self.dfast, self.lag, self.dlag
        ax.errorbar(lag, fast, xerr=dlag, yerr=dfast, color='white')

        ax.set_title('Likelihood') 
           
        return plt.show()

#
# class Group:
#
#     def __init__(self, listM, **kwargs):
#         """
#         A collection of measurements (or data).
#
#         args
#         -----
#
#         listM      list, a list of Q objects (all on same grid)
#
#         keyword arguments
#         -----------------
#
#         weights     numpy array, user-defined weights (order and size must match list)
#
#
#         example
#         --------
#
#         >>> S = Group(listM)           # initiate a Stack of EigenM objects
#         >>> l1l2stack = S.stack()           # Lam1 / Lam2 stack
#         >>> WSstack = S.wolfe_silver()      # Stack using Wolfe and Silver method
#         >>> RHstack = S.restivo_helffrich() # Stack using Restivo and Helffrich method
#         """
#
#         # get info
#         self.listM = listM
#
#         # check all grids are the same
#
#
#         # check degs
#         self.degs = self.listM[0].degs
#
#         nlags = min([ m._lags.size for m in ms ])
#         self.lags = self.listM[0].lags
#
#         # check all have the same grids
#         for Q in self.listM:
#             if np.any(M.degs != self.degs):
#                 raise Exception('Inconsistent degs grid found, all surface must have same grids.')
#             if np.any(M.lags != self.lags):
#                 raise Exception('Inconsistent lags grid found, all surfaces must have same grids.')
#
#         # if weights provided check it is the right size
#         if 'weights' in kwargs:
#             if not isinstance(kwargs['weights'], np.ndarray):
#                 raise TypeError('weights should be a numpy array')
#             if len(self.listM) != kwargs['weights'].size:
#                 raise Exception('weights array size must equal listM length')
#             self.weights = kwargs['weights'] # save weights for easy recall
#
#     def stack(self, **kwargs):
#
#         # if 'snr' in kwargs and kwargs['snr'] == True:
#         stack = np.stack([a.pdf for a in self.listM])
#         wts = [a.snr() for a in self.listM]
        


    #
    # def stack(self,**kwargs):
    #     """
    #     Return stack of lam1 / lam2 with optional user-defined weights.
    #     """
    #     listS = [ (M.lam1-M.lam2) / M.lam2 for M in self.listM ]
    #     return _stack(listS, **kwargs)
    #
    # def stackpdf(self,**kwargs):
    #     """
    #     Return stack of lam1 / lam2 with optional user-defined weights.
    #     """
    #     listS = [ ( (M.lam1-M.lam2) / M.lam2) / np.sum( (M.lam1-M.lam2) / M.lam2) for M in self.listM ]
    #     return _stack(listS, **kwargs)
#
# # basic stacking routine
# def _stack(listSurfaces,**kwargs):
#     """
#     Average listSurfaces, optionally using weights.
#
#     **kargs
#
#     weights = numpy array of weights
#     ----------------------------------
#     """
#
#     # put everything into one giant numpy array
#     stack = np.stack(listSurfaces)
#
#     # perform the averaging
#     # note np.average takes an optional weights parameter so automatically handles weigths.
#     return np.average(stack, axis=0, **kwargs)


