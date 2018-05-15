#!/usr/bin/env python

import numpy as np
import obspy
from os.path import dirname, join
import seispy
import glob
from datetime import timedelta


def read_catalog(logpath, b_time, e_time, stla, stlo):
    eq_lst = []
    with open(logpath) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            date_now = obspy.UTCDateTime(''.join(line.split()[0:7]))
            evla = float(line.split()[7])
            evlo = float(line.split()[8])
            evdp = float(line.split()[9])
            mw = float(line.split()[10])
            dis = seispy.distaz(stla, stlo, evla, evlo).delta
            bazi = seispy.distaz(stla, stlo, evla, evlo).getBaz()
            if b_time <= date_now <= e_time:
                eq_lst.append([date_now, evla, evlo, evdp, mw, dis, bazi])
    return eq_lst


def match_eq(eq_lst, pathname, ref_comp='Z', suffix='SAC', offset=0, tolerance=210):
    ref_eqs = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))
    eqs = []
    for evt in eq_lst:
        for ref_sac in ref_eqs:
            tr = obspy.read(ref_sac)[0]
            if tr.stats.starttime - timedelta(seconds=offset+tolerance) <= evt[0] <= tr.stats.starttime + timedelta(seconds=-offset+tolerance):
                eqs.append(pathname, eq(tr.stats.starttime, evt))


class eq(object):
    def __init__(self, pathname, evttime, evtinfo, dateformat='%Y.%j.%H.%M.%S'):
        self.st = obspy.read(join(pathname, '*'+evttime.strftime(dateformat)+'.*.SAC'))
        self.otime = evtinfo[0]
        self.evla = evtinfo[1]
        self.evlo = evtinfo[2]
        self.evdp = evtinfo[3]
        self.mw = evtinfo[4]
        self.dis = evtinfo[5]
        self.bazi = evtinfo[6]


class rf(object):
    def __init__(self, pathname):
        self.st = obspy.read(pathname)

    def retrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=order)

    def 
