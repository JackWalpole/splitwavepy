.. _introduction:

****************************************************
Get started
****************************************************

First off, get yourself a functioning `Python <https://www.python.org/>`_ environment, I use `Anaconda <https://www.anaconda.com/download/#macos>`_.

Second, install **SplitWavePy** with the following terminal command.

``pip install splitwavepy``

.. tip::
	Ensure you have the latest release.
	
	``pip upgrade splitwavepy``

Now try the code.

.. nbplot:: 
	:include-source:

	import splitwavepy as sw
	m = sw.EigenM( fast=50, lag=1.9, delta=0.05, noise=0.04)
	m.plot()

Save your measurement to disk.

.. nbplot::
	:include-source:
	
	>>> m.save('temp.eigm')


Check your working directory for a file called ``temp.eigm``.  Is it there?  How big is it?  It should be less than 100K.  It's a backup of your measurement together with the input data.

.. hint::
	To check your file size in a Unix style operating system:
	``du -sh temp.eigm``
	
Without closing your python session, load the data from the disk and check it's the same as that on memory.

.. nbplot::
	:include-source:
	
	>>> n = sw.load('temp.eigm')
	>>> n == m
	... True
	>>> n is m
	... False

	
Try plotting ``n``.  Does it look the same as ``m``?
	
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




