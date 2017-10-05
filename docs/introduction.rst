.. _introduction:

****************************************************
Introduction
****************************************************

Shear wave splitting is becoming an increasingly popular method to study anisotropy in the Earth, splitting is being measured in increasingly complex scenarios.

SplitWavePy provides simple to use open source tools to measure shear wave splitting on seismic datasets.  It works with either 2 (shear plane) -- or 3 -- component data.  Under the hood it's all numpy, which provides a quick gateway to fully unleashed Python.

Get started
------------

Import code, create a synthetic, measure shear-wave spliting, and plot, in just three lines:

.. nbplot:: 
	:include-source:
	
	import splitwavepy as sw
	m = sw.EigenM( fast=50, lag=1.9, delta=0.1, noise=0.03)
	m.plot()

This measures splitting using the powerful eigenvalue method of `Silver and Chan (1991) <http://onlinelibrary.wiley.com/doi/10.1029/91JB00899/abstract>`_.  

Making a measurement is as easy as instantiating an *EigenM* object.  If no arguments are supplied then the code will automatically produce a synthetic.  Check out the tutorials to see how to use `real data`_.

Keeping things together
-------------------------

Each measurement can be saved and backed up.  All information stored in the *EigenM* objected is preserved including the input data, any corrections that were applied as part of the measurement, and the :math:`\lambda_1` and :math:`\lambda_2` surfaces.

Saving and reloading the data is as easy as:

.. nbplot::
	
	m.save('temp.eigm')
	n = sw.load('temp.eigm')
	n == m

.. warning::

   Saving will overwrite pre-existing files with the same name.

From the loaded we can have a look at the original input data.

.. nbplot::
	:include-source:
	
	n.data.plot(window=True)
	
Or compare the :math:`\lambda_1` and :math:`\lambda_2` surfaces.

.. .. nbplot::
	:include-source:



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


Bootstrap error estimation 
--------------------------
 (Sandvol and Hearn)


Error surface stacking
----------------------


Lambda1 / Lambda2 surface stacking
----------------------------------


Contribute!
-----------

Code collaborators and pilot users, get in touch!
Please consider contributing to the code on github.

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




