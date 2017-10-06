.. _tutorial:

****************************************************
Tutorial
****************************************************

Assuming you have the code setup it's time to see what it can do.  Fire up an interactive python session such as ``ipython`` and ``import splitwavepy as sw``.

.. nbplot::
	:include-source:
	
	import splitwavepy as sw

Synthetic data
---------------------

By default, if no data is provided, SplitWavePy will make some synthetic data.  Data is stored in a *Pair* object.

.. .. autoclass:: splitwavepy.core.pair.Pair

.. nbplot::
	:include-source:

	data = sw.Pair(fast=30, lag=1.4, noise=0.03, pol=-29, delta=0.1)
	data.plot()

You specify the *fast* direction, *lag* time, *noise*, and source *pol* -arisation as shown above.  Order is not important.

.. note::
    Don't forget to set the sample interval *delta* appropriately (defaults to 1), it determines how your *lag* time is interpreted. 

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
    Every time a lag operation is applied the traces are shortened.  This is fine so long as your traces extend far enough beyond your analysis *Window*.

Window picking
----------------------------
	
The window should be designed in such a way as to maximise the energy of a single pulse of shear energy relative to noise and other arrivals.

Interactive plotting is supported by ``plot(interactive=True)``.  Click and drag to pick a window, or mouse hover and use the ``a`` and ``f`` keys to open and close the window (sac, anyone?), use the arrow keys to fine tune the window.  When you're ready to measure, hit the space bar.



*Windows* are parameterised by two (or three) parameters:

- *centre* time centre of window,
- *width* time length of window,
- *tukey* (optional) fraction of window to cosine taper (from 0 to 1).

.. note::
    You don't *need* to set the *Window*, the code will make a guess for you (designed to be about right for the synthetic case), but it's not very clever, so it's best to set it yourself.


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
	
Silver and Chan (1991) eigenvalue method
-----------------------------------------

A powerful and popular method for measuring splitting is the eigenvalue method of `Silver and Chan (1991) <http://onlinelibrary.wiley.com/doi/10.1029/91JB00899/abstract>`_.



How to do it
``````````````

With a window selected we are almost ready to meausure shear wave splitting!  We can tell the algorithm which splitting operators to trial using the *degs* and *tlags* keywords.  The measurement is made by instantiating an *EigenM* object.

.. .. autoclass:: splitwavepy.EigenM

>>> m = 60 # default 
>>>
>>> measurement = sw.EigenM( data, tlags=(mint,maxt,n))

.. note::

	If *Window*, *tlags*, or *degs* are unspecified, guesses are made.  It is strongly advised that you set these manually and at the very least check that these parameters look reasonable!
	
Making a measurement is as easy as instantiating an *EigenM* object.  If no arguments are supplied then the code will automatically produce a synthetic.  Check out the tutorials to see how to use :ref:`real_data`.


How it works
``````````````

Error Estimation
-----------------

F--test
````````

Bootstrap
``````````


The Results
----------------

.. Keeping things together
.. -------------------------
..
.. Each measurement can be saved and backed up to disk.
..
.. Saving and reloading the data is as easy as:
..
.. .. nbplot::
..
.. 	>>> m.save('temp.eigm')
.. 	>>> n = sw.load('temp.eigm')
.. 	>>> n == m
.. 	... True
..
.. .. warning::
..
..    Saving will overwrite pre-existing files with the same name.
..
..
.. **All** information stored in an *EigenM* objected is preserved, this includes:
..
.. * the input data,
.. * any corrections that were applied as part of the measurement, and
.. * the :math:`\lambda_1` and :math:`\lambda_2` surfaces.
..
.. From the loaded object we can look at the original input data.
..
.. .. nbplot::
.. 	:include-source:
..
.. 	n.data.plot()
..
.. Or compare the :math:`\lambda_1` and :math:`\lambda_2` surfaces.
..
.. .. nbplot::
.. 	:include-source:
..
.. 	fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(20, 6))
.. 	n.plot( ax=ax[0], vals=n.lam1, title=r'$\lambda_1$', mode='surf')
.. 	n.plot( ax=ax[1], vals=n.lam2, title=r'$\lambda_2$', mode='surf', cmap='magma_r')
.. 	n.plot( ax=ax[2], mode='surf') # by default plots (lam1-lam2)/lam2


An "error surface", in the literature, is typically a :math:`\lambda_2` surface with the values normalised so that the value at the confidence level of 95% is equal to 1.

- Error surface from F--test (using the summation coefficients found by Walsh et al., 2014).

>>> m.plot(m.error)

- Lambda1 and Lambda2 surfaces

>>> m.plot(m.lam1/m.lam2)

- Tabular data is reported with the following information:

+------------+------------+-----------+ 
| Header 1   | Header 2   | Header 3  | 
+============+============+===========+ 
| body row 1 | column 2   | column 3  | 
+------------+------------+-----------+ 

With optional geometry information appended to the right:

+------------+------------+-----------+ 
| Header 1   | Header 2   | Header 3  | 
+============+============+===========+ 
| body row 1 | column 2   | column 3  | 
+------------+------------+-----------+ 


Saving and loading
-------------------



.. _real_data:

Real data
---------

If you've got real data you need to get it into a `numpy <http://www.numpy.org/>`_ array.  `Obspy <https://github.com/obspy/obspy/wiki>`_ is extremely useful for that.  For the purposes of this tutorial, let's use obspy to download some data from the `<IRIS <https://www.iris.edu/hq/>`_ servers.


>>> import obspy

With real data it's worth doing a bit of pre-processing which at minimum will involve removing the mean from data, and might also involve bandpass filtering, interpolation, or rotating the components.  All of this is achievable in Obspy.

>>> # remove mean etc.


Once we're happy we can simply measure splitting by putting the data into a *Pair* and using the *EigenM* class as before

Transverse minimisation method
-------------------------------

Rotation correlation method
----------------------------

 
Null detection
--------------


Error surface stacking
----------------------


Self normalised SNR :math:`(\lambda_1 - \lambda_2)/\lambda_2` surface stacking
````````````````````````````````````````````````````````````````````````````````

If :math:`\lambda_1 = \text{signal} + \text{noise}` and :math:`\lambda_2 = \text{noise}`, then the signal to noise ratio, :math:`\text{SNR} = (\lambda_1 - \lambda_2)/\lambda_2`. 



Receiver correction
-------------------

Source correction
-----------------


Bootstrap correction error estimation
-------------------------------------






3--component data
--------------------







