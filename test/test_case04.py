import os
from seispy.rf import RF
from obspy import UTCDateTime
import pytest
import glob

def test_sub01():
    rf = RF()
    rf.para.datapath = './data'
    rf.para.phase = 'P'
    rf.para.stainfo.network = 'IC'
    rf.para.stainfo.station = 'BJT'
    rf.para.stainfo.channel = 'HH?'
    rf.para.magmin = 6.5
    rf.para.gauss = [1.0, 1.5]
    rf.para.use_remote_data = True
    rf.load_stainfo()

    rf.para.date_begin = UTCDateTime('20220101')
    rf.para.date_end = UTCDateTime('20220301')
    rf.search_eq()
    rf.match_eq()
    rf.save_raw_data()
    rf.detrend()
    rf.filter()
    rf.cal_phase()
    rf.rotate()
    rf.trim()

    rf.deconv()
    for ff in rf.para.gauss:
        rf.para.rfpath = './F{:.1f}/{}.{}'.format(ff, rf.para.stainfo.network, rf.para.stainfo.station)
        rf.saverf(ff)
    return rf


if __name__ == '__main__':
    test_sub01()