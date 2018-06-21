#!/usr/bin/env python

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

# save data
for tr in st:
    tr.write(tr.id + ".SAC", format="SAC")