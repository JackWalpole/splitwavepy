.. _transmin:

****************************************************
Transverse Minimisation Method
****************************************************

If the polarisation of the shear energy is known then the energy on the transverse component is indicative of splitting (Silver and Chan, 1998; Vinnik et al., 1989).  The method works by searching for the splitting parameters that when applied to the data removes the most energy from the transverse component.  In principle the method is very similar to the eigenvalue method of Silver and Chan 1991, except the polarisation of the shear energy us defined by the user, rather than being searched for by principal component analysis.  By including the polarisation as a priori information the measurement is more robust to noise, however, the method is vulnerable to error if the polarisation is specified incorrectly.

Continuing with the :ref:`realdata`, for which we obtained the following plot using the eigenvalue method.

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
	
To do the transverse minimisation method we use the ``TransM`` class with the polarisation specified.  In the case of the SKS wave the polarisation is equal to the backazimuth direction, this can be calculated using obspy.

.. nbplot::
	:include-source: True
	
	# use obspy to get backazimuth
	from obspy import geodetics
	dist, az, baz = geodetics.base.gps2dist_azimuth(evlat,evlon,stlat,stlon)
	
	# make the measurement
	m = sw.TransM(realdata, pol=baz, lags=(2,))
	m.plot()
	
Notice that the transverse minimisation method returns a more focussed result.

Now checkout the :ref:`rotcorr`.