.. _extensions:

****************************************************
Tutorial
****************************************************

Synthetic data
---------------------

Let's get started with some synthetic data.
This is achieved by the *Pair* class.

.. autoclass:: splitwavepy.core.pair.Pair

>>> # Start by importing
>>> import splitwavepy as sw
>>> import matplotlib.pyplot as plt
>>>
>>> # now it's a one-liner (set some parameters using keywords)
>>> data = sw.Pair(fast=30, lag=1.4, delta=0.1, noise=0.03)

You can change the *noise* level using the *noise* keyword.
Change the *fast* direction and the *lag* time similarly.
Don't forget to set the sample interval *delta* appropriately (defaults to 1) as it determines how your *lag* time is interpreted.  If you want to change the source polarisation use the *pol* keyword.

The *Pair* instance has some useful methods, one of which is *plot()*:

>>> data.plot()

Add another 'layer' of splitting.

>>> data.split( -12, 1.3) 
>>> data.plot()

Sometimes we might want to do a *correction* and apply the *inverse* splitting operator.

>>> data.unsplit( -12, 1.3)
>>> data.plot()

.. note::
	Every time a lag operation is applied the traces are shortened.
	
When measuring splitting we need a *Window*, which is defined by two parameters:

- *offset* from centre,
- *width* of the window,
- *tukey* (optional, 0 to 1) fraction of window to cosine taper.

.. autoclass:: splitwavepy.core.window.Window

.. note::

	This brings me to a subtle point about SplitWavePy, it works by a *centrality* principle.  Lags maintain balance on the centre sample and so must always have an even number of samples (x trace shifts half *lag* to the left, y trace shifts half *lag* to the right).  To ensure a balanced centre point all *Window* objects must have an odd *width*.  This should affect how you pick a *Window*.  You want the shear energy  in the middle of the *Window*, narrow enough not to incorporate too much surrounding energy, and wide enough to allow the energy to *spread* with lagging.
	
>>> wind = data.getWindow()
>>> data.plot(window=True)

With a window selected we are almost ready to meausure shear wave splitting!  We can tell the algorithm which splitting operators to trial using the *degs* and *tlags* keywords.  The measurement is made by instantiating an *EigenM* object.

.. autoclass:: splitwavepy.EigenM

>>> import numpy as np
>>> search_degs = np.linspace()
>>> search_tlags = np.linspace()
>>> m = sw.EigenM( data, degs=search_degs, tlags=search_tlags)

.. note::

	If *Window*, *tlags*, or *degs* are unspecified a guess is made.
	

Real data
---------

An easy way to access real data is by downloading it via Obspy.  In principle you can use any data so long as you can get it into a numpy array.

>>> import obspy

With real data it's worth doing a bit of pre-processing which at minimum will involve removing the mean from data, and might also involve bandpass filtering, interpolation, or rotating the components.  All of this is achievable in Obspy.

>>> # remove mean etc.


Once we're happy we can simply measure splitting by putting the data into a *Pair* and using the *EigenM* class as before




F--test error estimation
------------------------

We follow the method of 


Receiver correction
-------------------

Source correction
-----------------


Bootstrap correction error estimation
-------------------------------------



Null detection
--------------


3--component data
--------------------







