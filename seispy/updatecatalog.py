from datetime import datetime
from os.path import dirname, join
import urllib.request as rq
import re
import argparse
import sys
import os


def fetch_cata(inlog=join(dirname(__file__), 'data', 'EventCMT.dat'), outlog=''):
    url = 'http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk'
    try:
        print('Connecting to http://www.ldeo.columbia.edu/.../qcmt.ndk')
        response = rq.urlopen(url)
    except:
        raise TimeoutError('Could not connect to http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk')

    html = str(response.read())
    find_re = re.compile(r'[A-Z]+\s\d+/\d+/\d+\s.+?\\n', re.DOTALL)
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

    print('Writing station info to ' + outlog)
    for info in find_re.findall(html):
        year = int(info.split()[1].split('/')[0])
        mon = int(info.split()[1].split('/')[1])
        day = int(info.split()[1].split('/')[2])
        hour = int(info.split()[2].split(':')[0])
        min = int(info.split()[2].split(':')[1])
        sec = float(info.split()[2].split(':')[2])
        sec = round(sec)
        if sec == 60:
            sec = 59
        lat = float(info.split()[3])
        lon = float(info.split()[4])
        dep = info.split()[5]
        mw = info.split()[7]
        evt_time = datetime(year, mon, day, hour, min, sec)
        if old_time_end < evt_time:
            fid_new.write('%d %d %d %s %d %d %d %6.2f %6.2f %s %s\n' % (
            year, mon, day, evt_time.strftime('%j'), hour, min, sec, lat, lon, dep, mw))
    fid_new.close()


def main():
    parser = argparse.ArgumentParser(description="Update CMT Catalog")
    parser.add_argument('-i', help='Input Catalog', dest='inlog',
                        type=str, default=join(dirname(__file__), 'data', 'EventCMT.dat'))
    parser.add_argument('-o', help='Onput Catalog', dest='outlog', type=str, default='')
    arg = parser.parse_args()
    fetch_cata(inlog=arg.inlog, outlog=arg.outlog)


if __name__ == '__main__':
    pass