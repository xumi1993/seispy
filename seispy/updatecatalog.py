from datetime import datetime
from os.path import dirname, join
import urllib.request as rq
import re
import argparse
import sys
import os


def convertinfo(info):
    if info[5] == '60.0':
        sec = 59
    else:
        sec = float(info[5])
    return int(info[0]), int(info[1]), int(info[2]), int(info[3]), int(info[4]), sec, float(info[6]), float(info[7]), \
           float(info[8]), float(info[9])


def ndkparse(ndk_str):
    find_re = re.compile(
        r'[A-Z]+\s(\d+)/(\d+)/(\d+)\s+(\d+):(\d+):(\d+.\d)\s+(.+?)\s+(.+?)\s+(.+?)\s+.+?\s+(.+?)\s+.+?\n', re.DOTALL)
    ndk_lst = [convertinfo(info) for info in find_re.findall(ndk_str)]
    return ndk_lst


def fetch_cata(inlog=join(dirname(__file__), 'data', 'EventCMT.dat'), outlog=''):
    url = 'http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk'
    try:
        print('Connecting to http://www.ldeo.columbia.edu/.../qcmt.ndk')
        response = rq.urlopen(url)
    except Exception as e:
        raise TimeoutError('Could not connect to http://www.ldeo.columbia.edu/~gcmt/projects/CMT/'
                           'catalog/NEW_QUICK/qcmt.ndk\n{}'.format(e))

    html = response.read().decode('utf-8')
    with open(inlog, 'r') as fid_old:
        all_old_log = fid_old.readlines()
    old_log = all_old_log[-1]
    old_time_end = datetime(int(old_log.split()[0]), int(old_log.split()[1]), int(old_log.split()[2]),
                            int(old_log.split()[4]), int(old_log.split()[5]), int(old_log.split()[6]))

    if outlog == '':
        fid_new = open(inlog, 'w+')
        fid_new.writelines(all_old_log)
    else:
        fid_new = open(outlog, 'a+')

    print('Writing event info')
    for year, mon, day, hour, min, sec, lat, lon, dep, mw in ndkparse(html):
        evt_time = datetime(year, mon, day, hour, min, int(sec))
        if old_time_end < evt_time:
            fid_new.write('%d %d %d %s %d %d %d %6.2f %6.2f %s %s\n' % (
            year, mon, day, evt_time.strftime('%j'), hour, min, int(sec), lat, lon, dep, mw))
    fid_new.close()


def main():
    parser = argparse.ArgumentParser(description="Update CMT Catalog")
    parser.add_argument('-i', help='Input Catalog', dest='inlog', metavar='input_catalog',
                        type=str, default=join(dirname(__file__), 'data', 'EventCMT.dat'))
    parser.add_argument('-o', help='Onput Catalog', dest='outlog', metavar='output_catalog', 
                        type=str, default=join(dirname(__file__), 'data', 'EventCMT.dat'))
    # parser.add_argument('-u', help='url of ndk', metavar='ndl_url', type=str)
    arg = parser.parse_args()
    fetch_cata(inlog=arg.inlog, outlog=arg.outlog)


def ndk2dat():
    parser = argparse.ArgumentParser(description="Convert ndk file to dat catalog")
    parser.add_argument('-i', help='Input Catalog', dest='inlog', metavar='input_catalog',
                        type=str)
    parser.add_argument('-o', help='Onput Catalog', dest='outlog', metavar='output_catalog', 
                        type=str, default='EventCMT.dat')
    arg = parser.parse_args()
    with open(arg.inlog) as f:
        content = f.read()
    try:
        log_lst = ndkparse(content)
    except Exception as e:
        raise TypeError('Error type for ndk file\n{}'.format(e))
    with open(arg.outlog, 'w+') as f:
        for year, mon, day, hour, min, sec, lat, lon, dep, mw in log_lst:
            evt_time = datetime(year, mon, day, hour, min, int(sec))
            f.write('%d %d %d %s %d %d %d %6.2f %6.2f %s %s\n' % (
            year, mon, day, evt_time.strftime('%j'), hour, min, int(sec), lat, lon, dep, mw))


if __name__ == '__main__':
    ndk2dat()
