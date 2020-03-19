import obspy
from seispy.rf import *
from seispy.updatecatalog import fetch_cata


def subTC1():
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    rfproj = RF()
    rfproj.date_begin = date_begin
    rfproj.date_end = date_end
    rfproj.search_eq()
    print(rfproj.eq_lst)


def subTC2():
    fetch_cata()


if __name__ == '__main__':
    subTC1()
    subTC2()
