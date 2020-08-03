import obspy
from seispy.rf import *
from seispy.updatecatalog import fetch_cata
from obspy import UTCDateTime


def setpar(rf):
    rf.para.datapath = '9F33'
    rf.para.rfpath = 'RF_9F33'
    rf.para.suffix = 'sac'
    rf.para.ref_comp = '101'
    rf.para.magmin = 5
    rf.para.time_before = 10
    rf.para.time_after = 90
    rf.para.offset = 0
    rf.para.tolerance = 60
    rf.para.date_begin = UTCDateTime('20110301')
    rf.para.date_end = UTCDateTime('20110307')
    rf.para.criterion = 'crust'


def subTC1():
    rf = RF()
    setpar(rf)
    rf.load_stainfo()
    rf.search_eq()
    rf.match_eq()
    rf.detrend()
    rf.filter()
    rf.cal_phase()
    rf.trim()
    rf.rotate()
    rf.deconv()
    rf.saverf()
    print(rf.eq_lst)


def subTC2():
    fetch_cata()


if __name__ == '__main__':
    subTC1()
    subTC2()
