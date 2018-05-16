#!/usr/bin/env python

import numpy as np
import obspy
import re
from os.path import dirname, join
import seispy
import glob
from datetime import timedelta
import pandas as pd
from obspy.taup import import TauPyModel


def datestr2regex(datestr):
    patten = datestr.replace('%Y', r'\d{4}')
    patten = patten.replace('%m', r'\d{2}')
    patten = patten.replace('%d', r'\d{2}')
    patten = patten.replace('%j', r'\d{2}')
    patten = patten.replace('%H', r'\d{2}')
    patten = patten.replace('%M', r'\d{2}')
    patten = patten.replace('%S', r'\d{2}')
    return patten


def read_catalog(logpath, b_time, e_time, stla, stlo):
    col = ['date', 'evla', 'evlo', 'evdp', 'mag', 'dis', 'bazi']
    eq_lst = pd.DataFrame(columns=col)
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
                eq_lst = eq_lst.append([[date_now, evla, evlo, evdp, mw, dis, bazi]], columns=col, ignore_index=True)
    return eq_lst


def load_station_info(pathname, ref_comp, suffix):
    ex_sac = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))[0]
    ex_tr = obspy.read(ex_sac)[0]
    return ex_tr.stats.network, ex_tr.stats.station, ex_tr.stats.sac.stla, ex_tr.stats.sac.stlo


def match_eq(eq_lst, pathname, ref_comp='Z', suffix='SAC', offset=0, tolerance=210, dateformat='%Y.%j.%H.%M.%S'):
    pattern = datestr2regex(dateformat)
    ref_eqs = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))
    eq_match = pd.DataFrame(column=eq_lst.columns)
    for i, evt in eq_lst:
        for ref_sac in ref_eqs:
            datestr = re.findall(pattern, ref_sac)[0]
            tr = obspy.read(ref_sac)[0]
            if tr.stats.starttime - timedelta(seconds=offset+tolerance) <= evt['date'] <= tr.stats.starttime + timedelta(seconds=-offset+tolerance):
                eqs.append(eq(pathname, datestr, evt))
                eq_match = eq_match.append(evt, ignore_index=True)


class eq(object):
    def __init__(self, pathname, datestr, evtinfo):
        self.st = obspy.read(join(pathname, '*'+datestr+'.*.SAC'))
        self.otime = evtinfo[0]
        self.evla = evtinfo[1]
        self.evlo = evtinfo[2]
        self.evdp = evtinfo[3]
        self.mw = evtinfo[4]
        self.dis = evtinfo[5]
        self.bazi = evtinfo[6]

    def raypara(model):
        """
        model: TauPyModel form obspy
        """



class para():
    def __init__(self):
        self.datapath = ''
        self.RFpath = ''
        self.offset = 0
        self.tolerance = 210
        self.dateformat = '%Y.%j.%H.%M.%S'
        self.date_begin = obspy.UTCDateTime('19700101')
        self.date_end = obspy.UTCDateTime.now()


class stainfo():
    def __init__(self):
        self.network = ''
        self.station = ''
        self.stla = 0.
        self.stlo = 0.

    def load_stainfo(self, pathname, ref_comp, suffix):
        (self.network, self.station, self.stla, self.stlo) = load_station_info(pathname, ref_comp, suffix)


class rf(object):
    def __init__(self, pathname, date_begin, date_end, ref_comp='BHZ', suffix='SAC', offset=0, tolerance=210, dateformat='%Y.%j.%H.%M.%S'):
        self.para = para()
        self.stainfo = stainfo()
        self.para.datapath = pathname
        self.para.offset = offset
        self.para.tolerance = tolerance
        self.para.dateformat = dateformat
        self.para.date_begin = date_begin
        self.para.date_end = date_end
        self.ref_comp = ref_comp
        self.suffix = suffix
        self.eq_lst = pd.DataFrame()
        
        self.stainfo.load_stainfo(self.pathname, ref_comp, suffix)
        
        self.model = TauPyModel('iasp91')

    def search_eq(self, logpath):
        self.eq_lst = read_catalog(logpath, self.date_begin, self.date_end, self.stainfo.stla, self.stainfo.stlo)

    def retrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=order)

    def 
