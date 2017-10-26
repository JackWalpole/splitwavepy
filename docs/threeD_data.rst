.. _3Cdata:

****************************************************
3--Component Data
****************************************************

SplitWavePy supports working with 3--component data and offers tools for rotating the data into a user-defined reference frame.  This is implemented by providing a *ray* direction, which ideally should be normal to the shear plane.

3--component data is easily incorporated into the eigenvalue method, offering additional information in a third eigenvalue (e.g. Walsh et al., 2014).  In principle, if a shear wave has been fully corrected for splitting, this should maximise the first eigenvalue and minimise the second and third eigenvalues.

3--component data is handled by the Trio class.  This works in much the same way as the Pair class.  If no data are provided as arguments, the class will produce synthetic data with parameters set by keywords.

.. nbplot::
	:include-source:
	
	import splitwavepy as sw
	
	a = sw.Trio( split=( -50, 2.4), delta=0.025, noise=0.01, ray=(20,30))
	a.plot()

The main complication in dealing with 3-components is choosing the axis about which to perform the rotations and lags that define the splitting operator.  As a rule, the code will use the axis defined by the ray.  A trio object must, therefore, contain a ray.  The ray is set using the ``ray=(azi,inc)`` keyword which expects a tuple containing the azimuth and incidence angle, by default ``ray=(0,0)``.  Note, that the ray can be any angle you like, for example you might want to have it pointing in the wavefront normal direction, which could be different from the geometric ray direction.  The ray is plotted as a red line on the particle motion plots.

You can rotate your data to the ray as follows.

.. nbplot::
	:include-source:
	
	a.rotate2ray()
	a.plot()
  
To measure splitting on a Trio use the Eig3dM class.  You don't need to rotate to the ray before doing this.  The code will internally ensure the data is rotated into the ray frame before grid searching for the splitting parameters.

.. nbplot::
	:include-source:
	
	b = sw.Eig3dM(a,lags=(3,))
	b.plot()
	
.. note::

	In 3-D *all* fast directions are defined in the **ray frame**.  This includes the directions you assign to synthetics, receiver- and source-corrections, and the fast direction you measure.  Directions are measured relative to "up", which is actually the projection of the true up direction in the plane normal to the ray, i.e., the *SV* direction.
	
.. note::

		When using the trio you must ensure your data x, y, z correspond to north, east, up directions.
	
