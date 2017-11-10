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

By default, if no data are provided, SplitWavePy will make some synthetic data.  Data is stored in a *Pair* object.
I will use synthetic data to demonstrate the basic features of the code.

.. .. autoclass:: splitwavepy.core.pair.Pair

.. nbplot::
	:include-source:

	data = sw.Pair( split=( 30, 1.4), noise=0.03, pol=-15.8, delta=0.05)
	data.plot()

You can simulate splitting by specifying the fast direction and lag time using the ``split = (fast, lag)`` argument, similarly you can specify the *noise* level, and source *pol* -arisation as shown in the example above.  Order is not important.

.. note::
    Don't forget to set the sample interval *delta* appropriately (defaults to 1), it determines how your *lag* time is interpreted. 

Adding and removing splitting
------------------------------

The data are already split, but that doesn't mean they can't be split again!   You can add splitting using the *split()* method.

.. nbplot::
	:include-source:
	
	data.split( -12, 1.3) # fast, lag 
	data.plot()
	
.. note::
	This method will also split the noise, which is probable undesirable.  A better way to make synthetic with multiple splitting events is to provide the ``split`` keyword with a list of parameters, e.g. ``split = [( 30, 1.4), ( -12, 1.3)]``, these will be applied in the order provided.

Sometimes we might want to do a *correction* and apply the *inverse* splitting operator.  This is supported using the *unsplit()* method.  If we correct the above data for the latter layer of anisotropy we should return the waveforms to their former state.

.. nbplot::
	:include-source:

	data.unsplit( -12, 1.3)
	data.plot()

.. note::  
	Every time a lag operation is applied the traces are shortened.  
	This is fine so long as your traces extend far enough beyond your window.  

.. _window:

Setting the window
----------------------------
	
The window should be designed in such a way as to maximise the energy of a single pulse of shear energy relative to noise and other arrivals.

Set the window using the ``set_window(start,end)`` method.

.. nbplot::
	:include-source:

	data.set_window( 15, 32) # start, end 
	data.plot()

.. note::
    By default the window will be centred on the middle of the trace with width 1/3 of the trace length, which is likely to be inappropriate, so make sure to set your window sensibly.
	 
.. tip::
	Interactive window picking is supported by ``plot(pick=True)``.  Left click to pick the window start, right click to pick the window end, hit the space bar to save and close.  If you don't want to save your window, simply quit the figure using one of the built in matplotlib methods (e.g. click on the cross in the top left corner or hit the ``q`` key).pwd

.. .. nbplot::
	:include-source:

.. .. tip::
.. 	If the interactive plotting is not working you might need to add ``backend : TkAgg`` as a line
.. 	in your ``~/.matplotlib/matplotlibrc`` file.

	
.. .. note::
..
.. 	This brings me to a subtle but fundamental point about SplitWavePy, it works by a *centrality* principle.  Every lag operation involves a shift in the data, and must maintain balance on the centre sample.  Therefore every shift must always be an even number of samples (x trace shifts half *lag* to the left, y trace shifts half *lag* to the right).  To ensure a balanced centre point all *Window* objects must have an odd *width*.  This should affect how you pick a *Window*.  You want the shear energy  in the middle of the *Window*, narrow enough to avoid surrounding energy, and wide enough to capture relevant energy with a bit extra for 'spreading room'.
	

	
Silver and Chan (1991) eigenvalue method
-----------------------------------------

A powerful and popular method for measuring splitting is the eigenvalue method of `Silver and Chan (1991) <http://onlinelibrary.wiley.com/doi/10.1029/91JB00899/abstract>`_.  It uses a grid search to find the inverse splitting parameters that best linearise the particle motion.  Linearisation is assessed by principal component analysis at each search node, taking the eigenvalues of the covariance matrix, where linearity maximises :math:`\lambda_1` and minimises :math:`\lambda_2`.  The code uses the ratio :math:`\lambda_1/\lambda_2` to find the best node (which is more stable than using only :math:`\lambda_1` or :math:`\lambda_2` as it accounts for the possibility that energy might be lost by sliding out of the window).

To use this method on your data.

.. nbplot::
	:include-source:
	
	measure = sw.EigenM(data)
	measure.plot()

.. Changing the surface display
.. ``````````````````````````````
..
.. It is quite common to plot the :math:`\lambda_2` surface.  To do this you use the keyword ``vals = measure.lam2``, in fact any combination of :math:`\lambda_1` and :math:`\lambda_2` can be plotted in this way, for example you might have noticed that by default the code plots :math:`(\lambda_1-\lambda_2)/\lambda_2`.  Additionally, the title can be changed ``title = r'$\lambda_2$'``, a marker can be added at the max :math:`\lambda_1/\lambda_2` location ``marker=True``, and the 95\% confidence contour can be plotted ``conf95=True``.
..
.. .. This latter is the contour at the value of :math:`\lambda_2` that is 95% of the time, according to an F-test, which uses the noise level on the corrected residual trace (:math:\lambda_2 min) and a data determined estimate of the degrees of freedom (the code uses the coefficients of `Walsh et al., 2014 <http://onlinelibrary.wiley.com/doi/10.1002/jgrb.50386/full>`_).  To change the colour bar use ``cmap`` to a valid matplotlib colourmap.
..
.. .. nbplot::
.. 	:include-source:
..
.. 	measure.plot(vals=measure.lam2, title=r'$\lambda_2$', marker=True, conf95=True, cmap='viridis_r')
	
.. This example demonstrates that the ratio of :math:`\lambda_1` to :math:`\lambda_2` is much more sharply focused around the solution than a single eigenvalue :math:`\lambda_2` or :math:`\lambda_1` surface.
..
.. Furthermore, :math:`\lambda_1/\lambda_2` is dimensionless, and automatically normalised to the signal to noise ratio.  It's energy is strongly focused in the 95% confidence region, as would be expected for a proper error surface.  These properties suggest (to me), that surfaces of this type are much more appropriate for error surface stacking techniques than the (scaled) :math:`\lambda_2` surfaces that are normally used.

.. _setgrid:

Setting the lag time grid search
``````````````````````````````````
The code automatically sets the maximum lag time to be half the window length.  To set the max search time manually you use the ``lags`` keyword.  This accepts a tuple of length 1, 2, or 3, and will be interpreted differently depending on this length.  The rules are as follows: for a 1-tuple ``lags = (maxlag,)``, a 2-tuple ``lags = (maxlag, nlags)``, and finally a 3-tuple ``(minlag, maxlag, nlags)``.  Alternatively will accept a numpy array containing all nodes to search.

Setting the fast direction grid search
````````````````````````````````````````

The code automatically grid searches every 2 degrees along the fast direction axis.  That's ``degs = 90`` nodes in total (180/2).  You can change this number using the ``degs`` keyword and providing an integer.  Alternatively will accept a numpy array containing all nodes to search.


.. Tabulating the result
.. ----------------------
..
.. Oftentimes it is useful to reduce your measurement to tabular form.
.. This facilitates the analysis of a set of measurements in a spreadsheet type environment.
.. This is achievable by the ``report()`` method.  By default this will print a
.. tabular summary of your measurement to screen.
..
.. - By default tabular data is reported with the following information:
..
.. +------------+------------+-----------+
.. | Header 1   | Header 2   | Header 3  |
.. +============+============+===========+
.. | body row 1 | column 2   | column 3  |
.. +------------+------------+-----------+
..
..
.. .. note::
.. 	By reducing your measurement to tabular form you are losing valuable information.  Therefore do not rely on ``report()`` to backup your measurements.
.. 	Backup your measurements using ``save()``.


Saving and loading your measurements
-------------------------------------

To save your measurement to disk simply use the ``save(filename)`` method.
This will backup the input data complete with the :math:`\lambda_1` and :math:`\lambda_2` surfaces.

This can be recovered at a later time using ``splitwavepy.load(filename)``.

This feature is demonstrated here :ref:`introduction`.

Splitting corrections
----------------------

In the case where you have a good estimation of the splitting parameters beneath the receiver or the source it is possible to correct the waveforms and to measure the residual splitting.  The residual splitting can then be attributed to anisotropy elsewhere along the path.

Let's consider a simple 2-layer case.

.. nbplot::
	:include-source:
	
	# srcside and rceiver splitting parameters
	srcsplit = (  30, 1.3)
	rcvsplit = ( -44, 1.7)
	
	# Create synthetic
	a = sw.Pair( split=([ srcsplit, rcvsplit]), noise=0.03, delta=0.02)

	# standard measurement
	m = sw.EigenM(a, lags=(3,))
	m.plot()
	
The *apparent* splitting measured above is some non-linear combination of the 2-layers (non-linear because the order of splitting is important).

Receiver correction
``````````````````````
If we know the layer 2 contribution we can back this off and resolve the splitting in layer 1 using the ``rcvcorr=(fast, lag)`` keyword.
	
.. nbplot::
	:include-source:
	
	m = sw.EigenM(a, lags=(3,), rcvcorr=(-44,1.7))
	m.plot()
	
If it's worked we should have measured splitting parameters of :math:`\phi=30` and :math:`\delta t =1.3`.
	
Source correction
``````````````````

Alternatively, if we know the layer 1 contribution we can use
``srccorr=(fast, lag)`` to correct for the source side anisotropy.

.. nbplot::
	:include-source:
	
	m = sw.EigenM(a, lags=(3,), srccorr=(30,1.3))	
	m.plot()
	
If this has worked we should have measured splitting parameters of :math:`\phi=-44` and :math:`\delta t =1.7`.

If we apply both the source and receiver correction to the above synthetic example we should yield a *null* result (no splitting).

.. nbplot::
	:include-source:
	
	m = sw.EigenM(a, lags=(3,), rcvcorr=(-44,1.7), srccorr=(30,1.3))
	m.plot()

We do as can be seen by the concentration of energy at delay time 0.


.. Measurement stacking
.. ---------------------













