import pandas as pd
from obspy import UTCDateTime
import numpy as np
import argparse
from obspy.clients.fdsn import Client
from netCDF4 import Dataset
import sys


def _cat2df(cat):
    cols = ['date', 'evla', 'evlo', 'evdp', 'mag', 'magtype']
    data = [[evt.origins[0].time, evt.origins[0].latitude, evt.origins[0].longitude, evt.origins[0].depth*0.001,
             evt.magnitudes[0].mag, evt.magnitudes[0].magnitude_type] for evt in cat if evt.origins[0].depth is not None]
    return pd.DataFrame(data, columns=cols)


def wsfetch(server, starttime=None, endtime=None, minlatitude=None,
            maxlatitude=None, minlongitude=None, maxlongitude=None,
            latitude=None, longitude=None, minradius=None,
            maxradius=None, mindepth=None, maxdepth=None,
            minmagnitude=None, maxmagnitude=None, magnitudetype=None,
            includeallorigins=None, includeallmagnitudes=None,
            includearrivals=None, eventid=None, limit=None, offset=None,
            orderby='time-asc', catalog=None, contributor=None):
    if not isinstance(server, str):
        raise TypeError('server name should be \'str\' type')
    locs = locals()
    locs.pop('server')
    client = Client(server)
    cat = client.get_events(**locs)
    cat_df = _cat2df(cat)
    return cat_df


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


def lsnc():
    parser = argparse.ArgumentParser(description="List all fields of netCDF file")
    parser.add_argument('-k', help='Key name of fields', type=str, dest='key', default=None)
    parser.add_argument('ncfile', type=str, help='Path to netCDF file')
    arg = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    ncdata = Dataset(arg.ncfile)
    if arg.key is None:
        print(ncdata.variables)
    else:
        print(ncdata.variables[arg.key])


if __name__ == '__main__':
    ncfile = '/Users/xumj/Researches/Tibet_MTZ/models/3D2017-09Sv-depth.nc'
    ncdata = Dataset(ncfile)
    nc2npz(ncdata, minlat=22, maxlat=40, minlon=80, maxlon=105, maxdep=900)
