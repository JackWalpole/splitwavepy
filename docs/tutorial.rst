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

When measuring splitting we need a *Window*, which is defined by two parameters:
- *offset* from centre, and 
- *width* of the window, optionally a
- *tukey* parameter from 0 to 1 defines the fraction of the window to cosine taper.

.. autoclass:: splitwavepy.core.window.Window

.. note::

	This brings me to a subtle point about SplitWavePy, it works by a *centrality* principle.  Lags maintain balance on the centre sample and so must always have an even number of samples.  To ensure a balanced centre point all *Window* objects must have an odd *width*.  This should affect how you pick a *Window*.  You want the shear energy  in the middle of the *Window*, narrow enough not to incorporate too much surrounding energy, and wide enough to allow the energy to *spread* with lagging.

Real data
---------



Receiver correction
-------------------

Source correction
-----------------


Bootstrap correction error estimation
-------------------------------------



Null detection
--------------


3 -- component data
--------------------







