from obspy import UTCDateTime
import argparse
from seispy.io import Query
import sys
import pandas as pd

def download_catalog(fname, server='IRIS', format='seispy', **kwargs):
    query = Query(server=server)
    query.get_events(**kwargs)
    write_catalog(query, fname, format)


def write_catalog(query, fname, format):
    if format == 'seispy':
        with open(fname, 'w') as f:
            for i, row in query.events.iterrows():
                f.write('{} {:.3f} {:.3f} {:.1f} {:.1f}\n'.format(
                    row['date'].strftime('%Y %m %d %j %H %M %S'),
                    row['evla'], row['evlo'], row['evdp'], row['mag']
                ))
    else:
        query.events_raw.write(fname, format=format)


def read_catalog_file(fname):
    col = ['date', 'evla', 'evlo', 'evdp', 'mag']
    colcata = ['year', 'month', 'day', 'hour', 'minute', 'second']
    data_cata = pd.read_table(fname, header=None, sep=' |\t', engine='python')
    date_cols = data_cata.loc[:, [0,1,2,4,5,6]]
    date_cols.columns = colcata
    dates = pd.to_datetime(date_cols)
    eq_lst = pd.concat([dates, data_cata.loc[:, 7:]], axis=1)
    eq_lst.columns = col
    return eq_lst


def main():
    parser = argparse.ArgumentParser(description="Download Catalog to local file")
    parser.add_argument('fname', help='File name of output catalog')
    parser.add_argument('-f', help='The file format to use (e.g. "QUAKEML"). '
                        'See https://docs.obspy.org/master/packages/autogen/obspy.core.event.Catalog.write.html#supported-formats'
                        ' for supported formats. In addition, the \'seispy\' is supported for writing to a text file',
                         metavar='format', default='seispy')
    parser.add_argument('-b', help='Start time following obspy.UTCDateTime format', default=None, metavar='starttime')
    parser.add_argument('-e', help='End time following obspy.UTCDateTime format', default=None, metavar='endtime')
    parser.add_argument('-d', help='Radial geographic constraints with <lat>/<lon>/<minradius>/<maxradius>',
                        metavar='<lat>/<lon>/<minradius>/<maxradius>', default=None)
    parser.add_argument('-m', help='Magnitude <minmagnitude>[/<maxmagnitude>]',
                        metavar='<minmagnitude>[/<maxmagnitude>]', default=None)
    parser.add_argument('-p', help='Focal depth <mindepth>[/<maxdepth>]', metavar='<mindepth>[/<maxdepth>]', default=None)
    parser.add_argument('-s', help='Service provider, defaults to IRIS', metavar='server', default='IRIS')
    parser.add_argument('-c', help='Catalog name, defaults to GCMT', metavar='catalog', default='GCMT')
    arg = parser.parse_args()
    args = {}
    if arg.c is not None:
        args['catalog'] = arg.c
    if arg.p is not None:
        try:
            values = [float(value) for value in arg.p.split('/')]
        except:
            raise ValueError('Error format with focal depth')
        if len(values) == 1:
            args['mindepth'] = values[0]
        elif len(values) == 2:
            args['mindepth'] = values[0]
            args['maxdepth'] = values[1]
        else:
            raise ValueError('Error format with focal depth')
    if arg.m is not None:
        try:
            values = [float(value) for value in arg.m.split('/')]
        except:
            raise ValueError('Error format with magnitude')
        if len(values) == 1:
            args['minmagnitude'] = values[0]
        elif len(values) == 2:
            args['minmagnitude'] = values[0]
            args['maxmagnitude'] = values[1]
        else:
            raise ValueError('Error format with focal depth')
    if arg.d is not None:
        try:
            values = [float(value) for value in arg.d.split('/')]
        except:
            raise ValueError('Error format with Radial geographic constraints')
        if len(values) == 4:
            args['latitude'] = values[0]
            args['longitude'] = values[1]
            args['minradius'] = values[2]
            args['maxradius'] = values[3]
        else:
            raise ValueError('Error format with radial geographic constraints')
    if arg.b is not None:
        try:
            args['starttime'] = UTCDateTime(arg.b)
        except:
            raise ValueError('-b: Error format with time string')
    if arg.e is not None:
        try:
            args['endtime'] = UTCDateTime(arg.e)
        except:
            raise ValueError('-e: Error format with time string')
    if args == {}:
        parser.print_usage()
        sys.exit(1)
    download_catalog(arg.fname, server=arg.s, format=arg.f, **args)

