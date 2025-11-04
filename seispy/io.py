import pandas as pd
from obspy import UTCDateTime, Catalog
from obspy.clients.fdsn import Client


def _cat2df(cat):
    cols = ['date', 'evla', 'evlo', 'evdp', 'mag', 'magtype']
    data = [[evt.origins[0].time, evt.origins[0].latitude, evt.origins[0].longitude, evt.origins[0].depth*0.001,
             evt.magnitudes[0].mag, evt.magnitudes[0].magnitude_type] for evt in cat if evt.origins[0].depth is not None]
    return pd.DataFrame(data, columns=cols)


class Query():
    def __init__(self, server='IRIS', **kwargs):
        self.client = Client(server, **kwargs)

    def get_events(self, starttime=None,
                   endtime=UTCDateTime.now(), 
                   **kwargs):
        """Get events from IRIS

        :param starttime: Start time of events, defaults to None
        :type starttime: :class:`obspy.UTCDateTime`, optional
        :param endtime: End time of events, defaults to UTCDateTime.now()
        :type endtime: :class:`obspy.UTCDateTime`, optional
        :return: Events
        :rtype: :class:`obspy.Catalog`
        """
        if endtime > UTCDateTime.now():
            endtime = UTCDateTime.now()
        events = Catalog()
        if endtime-starttime < 365 * 86400:
            events += self.client.get_events(starttime=starttime,
                                            endtime=endtime,
                                            orderby='time-asc', **kwargs)
        else:
            chunk_length = 365 * 86400
            while starttime <= endtime:
                if endtime-starttime < chunk_length:
                    nowtime = endtime
                else:
                    nowtime=starttime + chunk_length
                try:
                    events += self.client.get_events(starttime=starttime,
                                                endtime=nowtime,
                                                orderby='time-asc', **kwargs)
                except:
                    starttime += chunk_length
                    continue
                if starttime + chunk_length > endtime:
                    chunk = endtime - starttime
                    if chunk <= 1:
                        break
                starttime += chunk_length

        self.events = _cat2df(events)
        self.events_raw = events

    def write_events(self, filename, format='QUAKEML'):
        """Write events to file

        :param filename: Output filename
        :type filename: str
        :param format: Output format, defaults to 'QUAKEML'. See the ObsPy documentation for supported formats. 
                       The 'TXT' format is supported to write a space-separated text file.
        :type format: str, optional
        """
        if format.upper() == 'TXT':
            self.events.to_csv(filename, index=False, sep=' ', header=False)
        else:
            self.events_raw.write(filename, format=format)

    def get_stations(self, includerestricted=False, **kwargs):
        self.stations = self.client.get_stations(includerestricted=includerestricted, **kwargs)


