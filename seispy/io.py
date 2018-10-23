import obspy
import pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client


def _cat2df(cat):
    cols = ['date', 'evla', 'evlo', 'evdp', 'mag', 'magtype']
    data = [[evt.origins[0].time, evt.origins[0].latitude, evt.origins[0].longitude, evt.origins[0].depth*0.001,
             evt.magnitudes[0].mag, evt.magnitudes[0].magnitude_type] for evt in cat]
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


if __name__ == '__main__':
    wsfetch('IRIS', starttime=UTCDateTime(2018,10,20))