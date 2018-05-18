.. _realdata:

****************************************************
Real data example
****************************************************


If you've got real data you need to get it into a `numpy <http://www.numpy.org/>`_ array.    For the purposes of this tutorial, let's use `Obspy <https://github.com/obspy/obspy/wiki>`_  to download some data from the `IRIS <https://www.iris.edu/hq/>`_ servers.

Preparing data
````````````````````

With real data it's worth doing a bit of pre-processing which at minimum will involve removing the mean from data, and might also involve bandpass filtering, interpolation, and/or rotating the components.  It is also necessary to pick the shear arrival of interest.  All of this is achievable in Obspy.

.. nbplot::
	:include-source:

	from obspy import read
	from obspy.clients.fdsn import Client
	from obspy import UTCDateTime

	client = Client("IRIS")
	t = UTCDateTime("2000-08-15T04:30:0.000")
	st = client.get_waveforms("IU", "CCM", "00", "BH?", t, t + 60 * 60,attach_response=True)

	# filter the data (this process removes the mean)
	st.filter("bandpass",freqmin=0.01,freqmax=0.5)

	st.plot()

When measuring splitting we need to have a specific shear wave arrival to target.  Let's use the obspy taup module to home in on the SKS arrival.

.. nbplot::
	:include-source:

	from obspy.taup import TauPyModel

	# from location and time, get event information
	lat=-31.56
	lon=179.74

	# server does not accept longitude greater than 180.
	cat = client.get_events(starttime=t-60,endtime=t+60,minlatitude=lat-1,
	                  maxlatitude=lat+1,minlongitude=lon-1,maxlongitude=180)
	evtime = cat.events[0].origins[0].time
	evdepth = cat.events[0].origins[0].depth/1000
	evlat = cat.events[0].origins[0].latitude
	evlon = cat.events[0].origins[0].longitude

	# station information
	inventory = client.get_stations(network="IU",station="CCM",starttime=t-60,endtime=t+60)
	stlat = inventory[0][0].latitude
	stlon = inventory[0][0].longitude
	
	# get backazimuth (radial polarisation at station)
	from obspy import geodetics
	dist, az, baz = geodetics.base.gps2dist_azimuth(evlat,evlon,stlat,stlon)

	# find arrival times
	model = TauPyModel('iasp91')
	arrivals = model.get_travel_times_geo(evdepth,evlat,evlon,stlat,stlon,phase_list=['SKS'])
	skstime = evtime + arrivals[0].time

	# trim around SKS
	st.trim( skstime-30, skstime+30)
	st.plot()
	
	
Measuring splitting
````````````````````

Now we have prepared data we are ready to measure splitting.

First read the shear plane components (horizontals in this case) into a *Data* object.

.. nbplot::

	import splitwavepy as sw

	# get data into Pair object and plot
	north = st[1].data
	east = st[0].data
	sample_interval = st[0].stats.delta
	realdata = sw.Data(north, east, delta=sample_interval)
	realdata.plot()
	
.. note::

	Order is important.  SplitWavePy expects the North component first.
	
By chance the window looks not bad.  If you want to change it see :ref:`window`.  For now let's press on with measuring the splitting.  I set the maximum delay time to search using the lags keyword ``lags=(2,)`` as explained in :ref:`setgrid`.

.. nbplot::
	:include-source:
	
	# Eigenvalue Metho
	realdata.EigenM(lags=(2,)).plot(name='Eigenvalue Method', ppm=True)
	realdata.TransM(lags=(2,), pol=baz).plot(name='Transverse Minimisation Method', ppm=True)
	realdata.XcorrM(lags=(2,)).plot(name='Cross-correlation Method', ppm=True)

SplitWavePy can also calculate the splitting intensity (Chevrot, 2000) from the data within the window.

.. nbplot::
	:include-source:
	
	print(realdata.splitting_intensity(pol=baz))


