.. _introduction:

****************************************************
Get started
****************************************************

First off, get yourself a functioning `Python <https://www.python.org/>`_ environment, I use `Anaconda <https://www.anaconda.com/download/#macos>`_.

Second, install **SplitWavePy** with the following terminal command.

``pip install splitwavepy``

.. sidebar::
	Ensure you have the latest release.
	
	``pip upgrade splitwavepy``

Now try the code.

.. nbplot:: 
	:include-source:

	import splitwavepy as sw
	m = sw.EigenM( fast=50, lag=1.9, delta=0.1, noise=0.03)
	m.plot()

Save your measurement to disk.

.. nbplot::
	:include-source:
	
	>>> m.save('temp.eigm')


Check your working directory for a file called ``temp.eigm``.  Is it there?  How big is it?
	
Without closing your kernel load the data from the disk and check it's the same as the original measurement on memory.

.. nbplot::
	:include-source:
	
	>>> n = sw.load('temp.eigm')
	>>> n == m
	... True
	>>> n is m
	... False

	
Try plotting ``n``.
	
If you've made it to here, great, you seem to have a working version of SplitWavePy.

Now check out the :ref:`tutorial`.



.. To do
.. -----
..
.. - Interactive *Window* picking
.. - Transverse energy minimization method.
.. - Rotation correlation method.
.. - *Q* calculation for null identification.
.. - Cluster *Window* analysis
.. - Frequency analysis
.. - Splitting intensity




