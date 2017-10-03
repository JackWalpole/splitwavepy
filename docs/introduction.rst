.. _extensions:

****************************************************
Introduction
****************************************************

SplitWavePy provides simple to use tools to measure shear wave splitting on seismic datasets.  It works with either 2 (shear plane) -- or 3 -- component data.  Under the hood it's all numpy, which provides a quick gateway to fully unleashed Python.

Easy to use
------------

Import code, create a synthetic, and plot it, in three lines:

>>> import splitwavepy as sw
>>> m = sw.EigenM( fast=50, lag=1.9, delta=0.05, noise=0.03)
>>> m.plot()

A Silver and Chan (1991) measurement
--------------------------------------

The code uses the eigenvalue method of `Silver and Chan (1991) <http://onlinelibrary.wiley.com/doi/10.1029/91JB00899/abstract>`_.  
It saves the following information:

- waveform data and analysis window

>>> m.data.plot(window=True)

- F-test error surface (using the summation coefficients found by Walsh et al., 2014).

>>> m.plot(m.error)

- lambda1 and lambda2 surfaces

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




To do list
------------

- Interactive *Window* picking
- Transverse energy minimization method.
- Rotation correlation method.
- *Q* calculation for null identification.
- Cluster *Windowow* analysis
- Frequency analysis
- Splitting intensity

Contribute!
-----------

Code collaborators and pilot users, get in touch!
Please consider contributing to the code on github.





