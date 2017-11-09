# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import core
from .window import Window

# class Traces:
#
#     def __init__(self, delta=None, *args, **kwargs):
#
#         # set delta
#         if delta is None: raise ValueError('delta must be set')
#         self.delta = delta
#
#         # set window
#         if 'window' in kwargs:
#             if isinstance(kwargs['window'], Window):
#                 self.__window = kwargs['window']
#                 return
#             else:
#                 raise TypeError('expecting a window')
#         # if no arguments provided
#         if len(args) == 0:
#             width = core.odd(self._nsamps() / 3)
#             self.__window = Window(width)
#             return
#         # if start/end given
#         if len(args) == 2:
#             start, end = args
#             self.__window = self._construct_window(start, end, **kwargs)
#             return
#         else:
#             raise Exception ('unexpected number of arguments')
#
#
#     @property
#     def delta(self):
#         return self.__delta
#
#     @delta.setter
#     def delta(self, delta):
#         self.__delta = delta
#
#     @property
#     def window(self):
#         return self.__window
#
#
#
#     # Window
#
#     def wbeg(self):
#         """
#         Window start time.
#         """
#         sbeg = self.__window.start(self._nsamps())
#         return sbeg * self.delta
#
#     def wend(self):
#         """
#         Window end time.
#         """
#         send = self.__window.end(self._nsamps())
#         return send * self.delta
#
#     def wwidth(self):
#         """
#         Window width.
#         """
#         return (self.__window.width-1) * self.delta
#
#     def _construct_window(self,start,end,**kwargs):
#         if start > end: raise ValueError('start is larger than end')
#         time_centre = (start + end)/2
#         time_width = end - start
#         tcs = core.time2samps(time_centre,self.delta)
#         offset = tcs - self._centresamp()
#         # convert time to nsamples -- must be odd (even plus 1 because x units of deltatime needs x+1 samples)
#         width = core.time2samps(time_width,self.delta,'even') + 1
#         return Window(width,offset,**kwargs)
            
            
    
class Measure():
    pass
    
#
#     def __init__(self, *args, **kwargs):
#
#         # inherit from Traces
#         Traces.__init__(self,*args,**kwargs)
#
#
#         self.
#
#
#
#
#     @property
#     def slags():
#         return self.__slags
#
#     @property
#     def degs():
#         return self.__degs
#
#     @property
#     def rcvcorr():
#         return self.__rcvcorr
#
#     @property
#     def srccorr():
#         return self.__srccorr
#
#     @slags.setter
#     def slags(self,**kwargs):
#         # LAGS
#         minlag = 0
#         maxlag = self.wwidth() / 4
#         nlags  = 40
#         if 'lags' not in kwargs:
#             lags = np.linspace( minlag, maxlag, nlags)))
#         else:
#             if isinstance(kwargs['lags'],np.ndarray):
#                 lags = kwargs['lags']
#             elif isinstance(kwargs['lags'],tuple):
#                 if len(kwargs['lags']) == 1:
#                     lags = np.linspace( minlag, kwargs['lags'][0], nlags)
#                 elif len(kwargs['lags']) == 2:
#                     lags = np.linspace( minlag,*kwargs['lags'])
#                 elif len(kwargs['lags']) == 3:
#                     lags = np.linspace( *kwargs['lags'])
#                 else:
#                     raise Exception('Can\'t parse lags keyword')
#             else:
#                 raise TypeError('lags keyword must be a tuple or numpy array')
#         # convert lags to samples (must be even)
#         self.__slags = np.unique( core.time2samps( lags, self.delta, mode='even')).astype(int)
#
#     @degs.setter
#     def degs(self,**kwargs):
#         # DEGS
#         mindeg = -90
#         maxdeg = 90
#         degs = 90
#         if 'degs' not in kwargs:
#             self.__degs = np.linspace( mindeg, maxdeg, degs, endpoint=False)
#         else:
#             if isinstance(kwargs['degs'],np.ndarray):
#                 self.__degs = kwargs['degs']
#             elif isinstance(kwargs['degs'],int):
#                 self.__degs = np.linspace( mindeg, maxdeg, kwargs['degs'], endpoint=False)
#             else:
#                 raise TypeError('degs must be an integer or numpy array')
#
#
#     @rcvcorr.setter
#     def rcvcorr(self,**kwargs):
#         # receiver correction
#         self.rcvcorr = None
#         if ('rcvcorr' in kwargs):
#             if not isinstance(kwargs['rcvcorr'],tuple): raise TypeError('rcvcorr must be tuple')
#             if len(kwargs['rcvcorr']) != 2: raise Exception('rcvcorr must be length 2')
#             # convert time shift to nsamples -- must be even
#             deg, lag = kwargs['rcvcorr']
#             samps = core.time2samps( lag,self.delta, 'even')
#             kwargs['rcvcorr'] = ( deg, samps)
#             self.rcvcorr = ( deg, samps * self.delta)
#
#     @srccorr.setter
#     def srccorr(self,**kwargs):
#         # source correction
#         self.srccorr = None
#         if ('srccorr' in kwargs):
#             if not isinstance(kwargs['srccorr'],tuple): raise TypeError('srccorr must be tuple')
#             if len(kwargs['srccorr']) != 2: raise Exception('srccorr must be length 2')
#             # convert time shift to nsamples -- must be even
#             deg, lag = kwargs['srccorr']
#             samps = core.time2samps( lag, self.delta, 'even')
#             kwargs['srccorr'] = ( deg, samps)
#             self.srccorr = ( deg, samps * self.delta)
#
#
#     def gridsearch(self,func):
#         """
#         Grid search for splitting parameters applied to data.
#         rcvcorr = receiver correction parameters in tuple (fast,lag)
#         srccorr = source correction parameters in tuple (fast,lag)
#         """
#
#         # grid of degs and lags to search over
#         degs, lags = np.meshgrid(self.degs,self.samplags)
#         shape = degs.shape
#         lam1 = np.zeros(shape)
#         lam2 = np.zeros(shape)
#
#         # avoid using "dots" in loops for performance
#         rotate = core.rotate
#         lag = core.lag
#         chop = core.chop
#         eigvalcov = core.eigvalcov
#
#         # ensure trace1 at zero angle
#         self.data.rotateto(0)
#         x, y = self.data.x, self.data.y
#
#         # pre-apply receiver correction
#         if 'rcvcorr' in kwargs:
#             x,y = core.unsplit(x,y,*kwargs['rcvcorr'])
#
#         # make function to do source correction (used in loop)
#         if 'srccorr' in kwargs:
#             srcphi, srclag = kwargs['srccorr']
#             def srccorr(x,y,ang):
#                 # unwind rotation
#                 x,y = rotate(x,y,srcphi-ang)
#                 # remove splitting
#                 x,y = lag(x,y,-srclag)
#                 return x,y
#         else:
#             def srccorr(x,y,ang):
#                 # no source correction so do nothing
#                 return x,y
#
#         for ii in np.arange(shape[1]):
#             tx, ty = rotate(x,y,degs[0,ii])
#             for jj in np.arange(shape[0]):
#                 # remove splitting so use inverse operator (negative lag)
#                 ux, uy = lag(tx,ty,-lags[jj,ii])
#                 # if requested -- post-apply source correction
#                 ux, uy = srccorr(ux,uy,degs[0,ii])
#                 # chop to analysis window
#                 ux, uy = chop(ux,uy,window=self.window)
#                 # measure eigenvalues of covariance matrix
#                 output[jj,ii] = eigvalcov(np.vstack((ux,uy)))
#
#         return output
        
