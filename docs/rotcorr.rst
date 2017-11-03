.. _rotcorr:

****************************************************
Rotation Correlation Method
****************************************************

The rotation correlation method searches for the rotation and lag which maximises the similarity of pulse shapes and aligns them in time.  The method was introduced by Fukao, 1984.

In my testing I have found it is very important to use the normalised correlation coefficient.

Again, continuing with the :ref:`realdata`, for which we have used the eigenvalue method and the transverse minimisation method plots:

.. nbplot::
	:include-source: False
	
	from obspy import read
	from obspy.clients.fdsn import Client
	from obspy import UTCDateTime
	from obspy.taup import TauPyModel

	client = Client("IRIS")
	t = UTCDateTime("2000-08-15T04:30:0.000")
	st = client.get_waveforms("IU", "CCM", "00", "BH?", t, t + 60 * 60,attach_response=True)

	# filter the data (this process removes the mean)
	st.filter("bandpass",freqmin=0.01,freqmax=0.5)

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

	# find arrival times
	model = TauPyModel('iasp91')
	arrivals = model.get_travel_times_geo(evdepth,evlat,evlon,stlat,stlon,phase_list=['SKS'])
	skstime = evtime + arrivals[0].time

	# trim around SKS
	st.trim( skstime-30, skstime+30)
	
	import splitwavepy as sw

	# get data into Pair object and plot
	north = st[1].data
	east = st[0].data
	sample_interval = st[0].stats.delta
	realdata = sw.Pair(north, east, delta=sample_interval)
	measure = sw.EigenM(realdata, lags=(2,))
	measure.plot()


.. nbplot::
	:include-source: False
	
	# use obspy to get backazimuth
	from obspy import geodetics
	dist, az, baz = geodetics.base.gps2dist_azimuth(evlat,evlon,stlat,stlon)
	
	# make the measurement
	m = sw.TransM(realdata, pol=baz, lags=(2,))
	m.plot()
	
	
To do the rotation correlation method we use the ``CrossM`` class as follows.

.. nbplot::
   :include-source: True
	
   m = sw.CrossM(realdata, lags=(2,))
   m.plot()
	

.. Wuestefeld and Bokelmann (2007) noticed that near null (no splitting) data produced systematically different results depending on the choice of method.
.. This was quantified in the parameter Q (Wuestefeld, et al., 2010), which gives `null` data negative values, and `good` data positive values.  I haven't got around to putting Q in the code yet, however

