.. _tutorial:

****************************************************
Tutorial
****************************************************

Synthetic data
---------------------

Let's get started with some synthetic data.
This is achieved by the *Pair* class.

.. .. autoclass:: splitwavepy.core.pair.Pair

.. nbplot::
	:include-source:
	
	# Before using the code we need to import it.
	import splitwavepy as sw
	
	# now it's a one-liner (set some parameters using keywords)
	data = sw.Pair(fast=30, lag=1.4, delta=0.1, noise=0.03)


You can change the *noise* level using the *noise* keyword, and the *fast* direction and *lag* time similarly.
Don't forget to set the sample interval *delta* appropriately (defaults to 1) as it determines how your *lag* time is interpreted.  If you want to change the source polarisation use the *pol* keyword.

The *Pair* instance has some useful methods, one of which is *plot()*:

.. nbplot::
	:include-source:
	
	data.plot()

The data are already split, but that doesn't mean they can't be split again! You can add more splitting using the *split()* method.

.. nbplot::
	:include-source:
	
	data.split( -12, 1.3) # fast, lag 
	data.plot()

Sometimes we might want to do a *correction* and apply the *inverse* splitting operator.  This is supported using the *unsplit()* method.  If we correct the above data for the latter layer of anisotropy we should return the waveforms to their former state.

.. nbplot::
	:include-source:

	data.unsplit( -12, 1.3)
	data.plot()

.. note::
	Every time a lag operation is applied the traces are shortened.  This is fine so long as your traces extend far beyond your analysis *Window*.
	
When measuring splitting we need a *Window*. This should be designed to capture the region of interest, only data within the *Window* will be subject to splitting analysis.

*Windows* are parameterised by two (or three) parameters:

- *centre* time centre of window,
- *width* time length of window,
- *tukey* (optional) fraction of window to cosine taper (from 0 to 1).

You don't *need* to set the *Window*, the code will make a guess for you (and the guess is designed to be about right for the synthetic case), but in general this guess could be wildly inappropriate, so it's best to keep a close eye and set it yourself.

.. .. autoclass:: splitwavepy.core.window.Window

A *Window* can be generate using the *getWindow(centre,width)* method, and you can see the *Window* on the data by using `plot(window=True)`.

.. nbplot::
	:include-source:

	wind = data.getWindow( 22, 15) # centre, width 
	data.plot(window=True)
	
	
.. .. note::

	This brings me to a subtle but fundamental point about SplitWavePy, it works by a *centrality* principle.  Every lag operation involves a shift in the data, and must maintain balance on the centre sample.  Therefore every shift must always be an even number of samples (x trace shifts half *lag* to the left, y trace shifts half *lag* to the right).  To ensure a balanced centre point all *Window* objects must have an odd *width*.  This should affect how you pick a *Window*.  You want the shear energy  in the middle of the *Window*, narrow enough to avoid surrounding energy, and wide enough to capture relevant energy with a bit extra for 'spreading room'.
	
.. .. nbplot::
	:include-source:



With a window selected we are almost ready to meausure shear wave splitting!  We can tell the algorithm which splitting operators to trial using the *degs* and *tlags* keywords.  The measurement is made by instantiating an *EigenM* object.

.. .. autoclass:: splitwavepy.EigenM

>>> import numpy as np
>>> search_degs = np.linspace()
>>> search_tlags = np.linspace()
>>> m = sw.EigenM( data, degs=search_degs, tlags=search_tlags)

.. note::

	If *Window*, *tlags*, or *degs* are unspecified, guesses are made.  It is strongly advised that you set these manually and at the very least check that these parameters look reasonable!
	

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







