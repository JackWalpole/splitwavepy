.. _realdata:

****************************************************
3-Component Data
****************************************************

Traditionally splitting is evaluated on 2-component data where those two components are oriented to capture energy in the shear plane.

In reality we often start with 3-component data and make a choice as to which 2-components to analyse.  For example it is common to use only the horizontal components when analysing steeply incident teleseismic shear arrivals.  SplitWavePy supports working with 3-component data and offers tools for rotating the data into a user-defined reference frame.  Furthermore, 3-component data is easily incorporated into the eigenvalue method, offering additional information in a third eigenvalue (e.g. Walsh, 2014).  In principle, if a shear wave has been fully corrected for splitting, this should maximise the first eigenvalue and minimise the second and third eigenvalues.

3-Component data is handled by the Trio class.  This works in much the same way as the Pair class.  If no data are provided as arguments, the class will produce synthetic data with parameters set by keywords.

.. nbplot::
	:include-source:
	
	import splitwavepy as sw
	
	a = sw.Trio( split=( -50, 2.4), delta=0.025, noise=0.01, ray=(20,30))
	a.plot()

The main complication in dealing with 3-components is choosing the axis about which to perform the rotations and lags that define the splitting operator.  As a rule, the code will use the axis defined by the ray.  A trio object must, therefore, contain a ray.  The ray is set using the ``ray=(azi,inc)`` keyword which expects a tuple containing the azimuth and incidence angle, by default ``ray=(0,0)``.  Note, this ``ray`` needn't conform to the traditional meaning of the word *ray*, for example in anisotropic media the *ray* may differ from the wavefront normal direction and you could use the ``ray`` keyword to rotate the data to the true wavefront.  The ray is plotted as a red line on the particle motion plots.
  
To measure splitting on a Trio use the Eig3dM class.

.. nbplot::
	:include-source:
	
	b = sw.Eig3dM(a,lags=(3,))
	b.plot()
	
.. note::

	In 3-D all fast directions are defined in the *ray frame*.

	When working in three dimensions the geographical reference frame (north, east, up) isn't particularly helpful. It doesn't deal very well with all possible ray angles.  It's better to scrap it and use the ray frame.  Rays can always be linked back to the geographical frame by their azimuth and incidence angle.  In the ray frame we rely on the "up" direction as a reference angle.  The "up" direction is not really the up direction, unless the ray is perfectly horizontal, it is actually the projection of the up direction in the plane normal to the ray. This is also called the *SV* direction.  What if the ray is perfectly vertical?  Then there is no "up" direction.  In this case, the ray azimuth is the default direction, so for example if the ray is set with azimuth of 0 and incidence of 0 ``ray=(0,0)`` then the reference angle is equivalent to north.
	
