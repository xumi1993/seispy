import obspy
from seispy.rf import *


def iris_fetch():
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    rfproj = rf()
    rfproj.date_begin = date_begin
    rfproj.date_end = date_end
    rfproj.search_eq()
    print(rfproj.eq_lst)


if __name__ == '__main__':
    iris_fetch()
