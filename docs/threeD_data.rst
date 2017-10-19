.. _realdata:

****************************************************
3-Component Data
****************************************************

Traditionally splitting is evaluated on 2-component data where those two components are oriented to capture energy in the shear plane.

In reality we often start with 3-component data and make a choice as to which 2-components to analyse.  For example it is common to use only the horizontal components when analysing steeply incident teleseismic shear arrivals.  SplitWavePy supports working with 3-component data and offers tools for rotating the data into a user-defined reference frame.  Furthermore, 3-component data is easily incorporated into the eigenvalue method, offering additional information in a third eigenvalue (e.g. Walsh, 2014).  In principle, if a shear wave has been fully corrected for splitting, this should maximise the first eigenvalue and minimise the second and third eigenvalues.

3-Component data is handled by the Trio class.  This works in much the same way as the Pair class, if no data are provided as arguments, the class will produce synthetic data.

.. nbplot::
	:include-source:
	
	import splitwavepy as sw
	
	a = sw.Trio( split=( -50, 2.4), delta=0.025, noise=0.01, ray=(20,30))
	a.plot()
	
Notice we have the ``ray=(azi,inc)`` keyword which expects a tuple containing the azimuth and incidence angle of the ray.

To measure splitting on a Trio use the Eig3dM class.

.. nbplot::
	:include-source:
	
	b = sw.Eig3dM(a,lags=(3,))
	b.plot()