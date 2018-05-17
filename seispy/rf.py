#!/usr/bin/env python

import numpy as np
import obspy
import re
from os.path import dirname, join, expanduser
import seispy
import glob
from datetime import timedelta
import pandas as pd
from obspy.taup import TauPyModel
import deepdish as dd

def datestr2regex(datestr):
    pattern = datestr.replace('%Y', r'\d{4}')
    pattern = pattern.replace('%m', r'\d{2}')
    pattern = pattern.replace('%d', r'\d{2}')
    pattern = pattern.replace('%j', r'\d{3}')
    pattern = pattern.replace('%H', r'\d{2}')
    pattern = pattern.replace('%M', r'\d{2}')
    pattern = pattern.replace('%S', r'\d{2}')
    return pattern


def read_catalog(logpath, b_time, e_time, stla, stlo, magmin=5.5, magmax=10, dismin=30, dismax=90):
    col = ['date', 'evla', 'evlo', 'evdp', 'mag', 'dis', 'bazi']
    eq_lst = pd.DataFrame(columns=col)
    with open(logpath) as f:
        lines = f.readlines()
        for line in lines:
            line_sp = line.strip().split()
            date_now = obspy.UTCDateTime.strptime('.'.join(line_sp[0:3])+'T'+'.'.join(line_sp[4:7]), '%Y.%m.%dT%H.%M.%S')
            evla = float(line_sp[7])
            evlo = float(line_sp[8])
            evdp = float(line_sp[9])
            mw = float(line_sp[10])
            dis = seispy.distaz(stla, stlo, evla, evlo).delta
            bazi = seispy.distaz(stla, stlo, evla, evlo).getBaz()
            if b_time <= date_now <= e_time and magmin <= mw <= magmax and dismin <= dis <= dismax:
                this_data = pd.DataFrame([[date_now, evla, evlo, evdp, mw, dis, bazi]], columns=col)
                eq_lst = eq_lst.append(this_data, ignore_index=True)
    return eq_lst


def load_station_info(pathname, ref_comp, suffix):
    ex_sac = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))[0]
    ex_tr = obspy.read(ex_sac)[0]
    return ex_tr.stats.network, ex_tr.stats.station, ex_tr.stats.sac.stla, ex_tr.stats.sac.stlo


def match_eq(eq_lst, pathname, ref_comp='Z', suffix='SAC', offset=0, tolerance=210, dateformat='%Y.%j.%H.%M.%S'):
    pattern = datestr2regex(dateformat)
    ref_eqs = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))
    sac_files = []
    for ref_sac in ref_eqs:
        datestr = re.findall(pattern, ref_sac)[0]
        tr = obspy.read(ref_sac)[0]
        sac_files.append([datestr, tr])
    new_col = ['data']
    eq_match = pd.DataFrame(columns=new_col)
    print('Matching SAC files')
    for i, evt in eq_lst.iterrows():
        tmp_datestr = []
        for datestr, tr in sac_files:
            if tr.stats.starttime - timedelta(seconds=offset+tolerance) <= evt['date'] <= tr.stats.starttime + timedelta(seconds=-offset+tolerance):
                tmp_datestr.append(datestr)
        if len(tmp_datestr) == 1:
            this_eq = eq(pathname, tmp_datestr[0])
            this_df = pd.DataFrame([[this_eq]], columns=new_col, index=[i])
            eq_match = eq_match.append(this_df)
    return pd.concat([eq_lst, eq_match], axis=1, join='inner')


class eq(object):
    def __init__(self, pathname, datestr):
        self.st = obspy.read(join(pathname, '*'+datestr+'.*.SAC'))
        self.rf = obspy.Stream()
    
    def detrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=order)


class para():
    def __init__(self):
        self.datapath = expanduser('~')
        self.RFpath = expanduser('~')
        self.offset = 0
        self.tolerance = 210
        self.dateformat = '%Y.%j.%H.%M.%S'
        self.date_begin = obspy.UTCDateTime('19760101')
        self.date_end = obspy.UTCDateTime.now()
        self.magmin = 5.5
        self.magmax = 10
        self.dismin = 30
        self.dismax = 90
        self.ref_comp = 'BHZ'
        self.suffix = 'SAC'

    def get_para(self):
        return self.__dict__


class stainfo():
    def __init__(self):
        self.network = ''
        self.station = ''
        self.stla = 0.
        self.stlo = 0.

    def get_stainfo(self):
        return self.__dict__

    def load_stainfo(self, pathname, ref_comp, suffix):
        (self.network, self.station, self.stla, self.stlo) = load_station_info(pathname, ref_comp, suffix)


class rf(object):
    def __init__(self):
        self.para = para()
        self.stainfo = stainfo()
        self.eq_lst = pd.DataFrame()
        self.eqs = pd.DataFrame()   
        self.model = TauPyModel('iasp91')
    
    @property
    def date_begin(self):
        return self.para.date_begin

    @date_begin.setter
    def date_begin(self, value):
        self.para.date_begin = value

    @property
    def date_end(self):
        return self.para.date_end

    @date_end.setter
    def date_end(self, value):
        self.para.date_end = value
    
    @property
    def datapath(self):
        return self.para.datapath

    @datapath.setter
    def datapath(self, value):
        self.para.datapath = value

    def load_stainfo(self):
        self.stainfo.load_stainfo(self.para.datapath, self.para.ref_comp, self.para.suffix)

    def search_eq(self, logpath):
        self.eq_lst = read_catalog(logpath, self.para.date_begin, self.para.date_end,
                                   self.stainfo.stla, self.stainfo.stlo,
                                   magmin=self.para.magmin, magmax=self.para.magmax,
                                   dismin=self.para.dismin, dismax=self.para.dismax)

    def match_eq(self):
        self.eqs = match_eq(self.eq_lst, self.para.datapath, ref_comp=self.para.ref_comp, suffix=self.para.suffix,
                            offset=self.para.offset, tolerance=self.para.tolerance,
                            dateformat=self.para.dateformat)

    def save(self, path=''):
        if path == '':
            path = '{0}.{1}.npz'.format(self.stainfo.network, self.stainfo.station)
        
        para = self.para.__dict__
        stainfo = self.stainfo.__dict__
        
        d = {'para':para, 'stainfo':stainfo, 'eq_lst':self.eq_lst, 'eqs':self.eqs}
        try:
            dd.io.save(path, d)
        except Exception as e:
            raise IOError(e)

    def load(self, path):
        try:
            fdd = dd.io.load(path)
        except Exception as e:
            raise IOError('Cannot read {0}'.format(path))
        
        try:
            self.para.__dict__.update(fdd['para'])
            self.stainfo.__dict__.update(fdd['stainfo'])
            self.eq_lst = fdd['eq_lst']
            self.eqs = fdd['eqs']
        except Exception as e:
            raise ValueError(e)

    def retrend(self):
        for i, row in self.eqs:
            row['data'].detrend()

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        for i, row in self.eqs:
            row['data'].filter(freqmin=freqmin, freqmax=freqmax, corners=order)


def InitRfProj():
    rfproj = rf()


if __name__ == '__main__':
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    logpath = '/Users/xumj/Codes/seispy/Scripts/EventCMT.dat'
    # logpath = '/home/xu_mijian/Codes/seispy/Scripts/EventCMT.dat'
    datapath = '/Users/xumj/Researches/test4seispy/data'
    # datapath = '/home/xu_mijian/xu_mijian/NJ2_SRF/data'
    # proj_file = '/home/xu_mijian/xu_mijian/NJ2_SRF/test.h5'
    proj_file = '/Users/xumj/Researches/test4seispy/test.h5'

    rfproj = rf()
    rfproj.load(proj_file)
    # rfproj.date_begin = date_begin
    # rfproj.date_end = date_end
    # rfproj.datapath = datapath
    # rfproj.load_stainfo()
    # rfproj.search_eq(logpath)
    # rfproj.match_eq()
    # rfproj.save(proj_file)



