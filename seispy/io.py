import pandas as pd
from obspy import UTCDateTime, Catalog
import numpy as np
from datetime import timedelta
from obspy.clients.fdsn import Client
# from netCDF4 import Dataset
import sys


def _cat2df(cat):
    cols = ['date', 'evla', 'evlo', 'evdp', 'mag', 'magtype']
    data = [[evt.origins[0].time, evt.origins[0].latitude, evt.origins[0].longitude, evt.origins[0].depth*0.001,
             evt.magnitudes[0].mag, evt.magnitudes[0].magnitude_type] for evt in cat if evt.origins[0].depth is not None]
    return pd.DataFrame(data, columns=cols)


class Query():
    def __init__(self, server='IRIS'):
        self.client = Client(server)

    def get_events(self, starttime=None,
                   endtime=UTCDateTime.now(), **kwargs):
        if endtime > UTCDateTime.now():
            endtime = UTCDateTime.now()
        chunk_length = 365 * 86400
        events = Catalog()
        while starttime <= endtime:
            events += self.client.get_events(starttime=starttime,
                                             endtime=starttime + chunk_length,
                                             **kwargs)
            if starttime + chunk_length > endtime:
                chunk = endtime - starttime
                if chunk <= 1:
                    break
            starttime += chunk_length
        self.events = _cat2df(events)

    def get_stations(self, includerestricted=False, **kwargs):
        self.stations = self.client.get_stations(includerestricted=includerestricted, **kwargs)


def nc2npz(ncdata, minlat=-90, maxlat=90, minlon=-180, maxlon=180, mindep=0, maxdep=6371, key='dvs'):
    lat = ncdata.variables['latitude'][:].data
    lon = ncdata.variables['longitude'][:].data
    dep = ncdata.variables['depth'][:].data
    data = ncdata.variables[key][:].data
    idx_lat = np.where((lat >= minlat) & (lat <= maxlat))[0]
    idx_lon = np.where((lon >= minlon) & (lon <= maxlon))[0]
    idx_dep = np.where((dep >= mindep) & (dep <= maxdep))[0]
    cut_data = data[idx_dep[0]:idx_dep[-1]+1, idx_lat[0]:idx_lat[-1]+1, idx_lon[0]:idx_lon[-1]+1]
    cut_lat = lat[idx_lat]
    cut_lon = lon[idx_lon]
    cut_dep = dep[idx_dep]
    # new_lat, new_dep, new_lon = np.meshgrid(cut_lat, cut_dep, cut_lon)
    return cut_data, cut_dep, cut_lat, cut_lon

