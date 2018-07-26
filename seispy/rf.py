#!/usr/bin/env python

import numpy as np
import obspy
import re
import io
from os.path import dirname, join, expanduser, exists
from seispy.para import para
import seispy
from seispy.eq import eq
from seispy.setuplog import setuplog
import glob
from datetime import timedelta, datetime
import pandas as pd
from pandas.errors import PerformanceWarning
from obspy.taup import TauPyModel
import deepdish as dd
import urllib.request as rq
import configparser


def datestr2regex(datestr):
    pattern = datestr.replace('%Y', r'\d{4}')
    pattern = pattern.replace('%m', r'\d{2}')
    pattern = pattern.replace('%d', r'\d{2}')
    pattern = pattern.replace('%j', r'\d{3}')
    pattern = pattern.replace('%H', r'\d{2}')
    pattern = pattern.replace('%M', r'\d{2}')
    pattern = pattern.replace('%S', r'\d{2}')
    return pattern


def get_events(b_time, e_time, stla, stlo, magmin=5.5, magmax=10, dismin=30, dismax=90):
    starttime = b_time.strftime('%Y-%m-%dT%H:%M:%S')
    endtime = e_time.strftime('%Y-%m-%dT%H:%M:%S')
    use_cols = ['Time', 'Latitude', 'Longitude', 'Depth', 'Magnitude']
    real_cols = ['date', 'evla', 'evlo', 'evdp', 'mag']
    dateparse = lambda x: obspy.UTCDateTime.strptime(x, '%Y-%m-%dT%H:%M:%SZ')
    url = 'http://service.iris.edu/fdsnws/event/1/query?&starttime={0}' \
          '&endtime={1}&lat={2}&lon={3}&minradius={4}&' \
          'maxradius={5}&minmag={6}&maxmag={7}&catalog=GCMT&' \
          'orderby=time-asc&format=geocsv'.format(starttime, endtime, stla, stlo, dismin, dismax, magmin, magmax)
    try:
        response = rq.urlopen(url)
    except Exception as e:
        raise ConnectionError('{0}'.format(e))
    evt_csv = io.StringIO(response.read().decode())
    df = pd.read_csv(evt_csv, sep='|', comment='#', parse_dates=[0], date_parser=dateparse, usecols=use_cols)
    # daz = seispy.distaz(stla, stlo, df['Latitude'].values, df['Longitude'].values)
    # df['dis'] = pd.Series(daz.delta, index=df.index)
    # df['bazi'] = pd.Series(daz.baz, index=df.index)
    col_dict = dict(zip(use_cols, real_cols))
    df.rename(columns=col_dict, inplace=True)
    return df


def read_catalog(logpath, b_time, e_time, stla, stlo, magmin=5.5, magmax=10, dismin=30, dismax=90):
    col = ['date', 'evla', 'evlo', 'evdp', 'mag']
    eq_lst = pd.DataFrame(columns=col)
    with open(logpath) as f:
        lines = f.readlines()
        for line in lines:
            line_sp = line.strip().split()
            date_now = obspy.UTCDateTime.strptime('.'.join(line_sp[0:3]) + 'T' + '.'.join(line_sp[4:7]),
                                                  '%Y.%m.%dT%H.%M.%S')
            evla = float(line_sp[7])
            evlo = float(line_sp[8])
            evdp = float(line_sp[9])
            mw = float(line_sp[10])
            dis = seispy.distaz(stla, stlo, evla, evlo).delta
            # bazi = seispy.distaz(stla, stlo, evla, evlo).getBaz()
            if b_time <= date_now <= e_time and magmin <= mw <= magmax and dismin <= dis <= dismax:
                this_data = pd.DataFrame([[date_now, evla, evlo, evdp, mw]], columns=col)
                eq_lst = eq_lst.append(this_data, ignore_index=True)
    return eq_lst


def load_station_info(pathname, ref_comp, suffix):
    try:
        ex_sac = glob.glob(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))[0]
    except Exception:
        raise FileNotFoundError('no such SAC file in {0}'.format(pathname))
    ex_tr = obspy.read(ex_sac)[0]
    return ex_tr.stats.network, ex_tr.stats.station, ex_tr.stats.sac.stla, ex_tr.stats.sac.stlo, ex_tr.stats.sac.stel


def match_eq(eq_lst, pathname, stla, stlo, ref_comp='Z', suffix='SAC', offset=None,
             tolerance=210, dateformat='%Y.%j.%H.%M.%S', switchEN=False, reverseE=False, reverseN=False):
    pattern = datestr2regex(dateformat)
    ref_eqs = glob.glob(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))
    sac_files = []
    for ref_sac in ref_eqs:
        datestr = re.findall(pattern, ref_sac)[0]
        try:
            tr = obspy.read(ref_sac)[0]
        except TypeError:
            continue
        if offset is None:
            this_offset = tr.stats.sac.o
        elif isinstance(offset, (int, float)):
            this_offset = -offset
        else:
            raise TypeError('offset should be int or float type')
        sac_files.append([datestr, tr.stats.starttime, this_offset])
    new_col = ['dis', 'bazi', 'data']
    eq_match = pd.DataFrame(columns=new_col)
    for datestr, b_time, offs in sac_files:
        date_range_begin = b_time + timedelta(seconds=offs - tolerance)
        date_range_end = b_time + timedelta(seconds=offs + tolerance)
        results = eq_lst[(eq_lst.date > date_range_begin) & (eq_lst.date < date_range_end)]
        if len(results) != 1:
            continue
        try:
            this_eq = eq(pathname, datestr, suffix, switchEN=switchEN, reverseE=reverseE, reverseN=reverseN)
        except Exception as e:
            continue
        this_eq.get_time_offset(results.iloc[0]['date'])
        daz = seispy.distaz(stla, stlo, results.iloc[0]['evla'], results.iloc[0]['evlo'])
        this_df = pd.DataFrame([[daz.delta, daz.baz, this_eq]], columns=new_col, index=results.index.values)
        eq_match = eq_match.append(this_df)
    ind = eq_match.index.drop_duplicates(False)
    eq_match = eq_match.loc[ind]
    '''
    for i, evt in eq_lst.iterrows():
        tmp_datestr = []
        for datestr, b_time, offs in sac_files:
            if b_time + timedelta(seconds=offs - tolerance) <= evt['date'] \
                    <= b_time + timedelta(seconds=offs + tolerance):
                tmp_datestr.append(datestr)
        if len(tmp_datestr) == 1:
            try:
                this_eq = eq(pathname, tmp_datestr[0], suffix)
            except:
                continue
            this_eq.get_time_offset(evt['date'])
            daz = seispy.distaz(stla, stlo, evt['evla'], evt['evlo'])
            this_df = pd.DataFrame([[daz.delta, daz.baz, this_eq]], columns=new_col, index=[i])
            eq_match = eq_match.append(this_df)
    '''
    return pd.concat([eq_lst, eq_match], axis=1, join='inner')


class stainfo():
    def __init__(self):
        self.network = ''
        self.station = ''
        self.stla = 0.
        self.stlo = 0.
        self.stel = 0.

    def get_stainfo(self):
        return self.__dict__

    def load_stainfo(self, pathname, ref_comp, suffix):
        (self.network, self.station, self.stla, self.stlo, self.stel) = load_station_info(pathname, ref_comp, suffix)


def CfgParser(cfg_file):
    cf = configparser.RawConfigParser()
    pa = para()
    try:
        cf.read(cfg_file)
    except Exception:
        raise FileNotFoundError('Cannot open configure file %s' % cfg_file)
    pa.datapath = cf.get('path', 'datapath')
    pa.rfpath = cf.get('path', 'rfpath')
    pa.imagepath = cf.get('path', 'imagepath')
    pa.catalogpath = cf.get('path', 'catalogpath')
    for key, value in cf.items('para'):
        if key == 'date_begin':
            pa.__dict__[key] = obspy.UTCDateTime(value)
        elif key == 'date_end':
            pa.__dict__[key] = obspy.UTCDateTime(value)
        elif key == 'offset':
            try:
                pa.__dict__[key] = float(value)
            except:
                pa.__dict__[key] = None
        elif key == 'only_r':
            pa.only_r = cf.getboolean('para', 'only_r')
        else:
            try:
                pa.__dict__[key] = float(value)
            except ValueError:
                pa.__dict__[key] = value
    return pa


class rf(object):
    def __init__(self, cfg=None, log=None):
        if cfg is None:
            self.para = para()
        elif isinstance(cfg, str):
            self.para = CfgParser(cfg)
        else:
            raise TypeError('cfg should be \'str\' not \'{0}\''.format(type(cfg)))
        if not isinstance(self.para, para):
            raise TypeError('Input value should be class seispy.rf.para')
        if log is None:
            self.logger = setuplog()
        else:
            self.logger = log
        self.eq_lst = pd.DataFrame()
        self.eqs = pd.DataFrame()
        self.model = TauPyModel('iasp91')
        self.stainfo = stainfo()

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

    def load_stainfo(self):
        try:
            self.logger.RFlog.info('Load station info from {0}'.format(self.para.datapath))
            self.stainfo.load_stainfo(self.para.datapath, self.para.ref_comp, self.para.suffix)
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise e

    def search_eq(self, local=False):
        if not local:
            try:
                self.logger.RFlog.info('Searching earthquakes from IRIS-WS')
                self.eq_lst = get_events(self.para.date_begin, self.para.date_end, self.stainfo.stla,
                                         self.stainfo.stlo, magmin=self.para.magmin, magmax=self.para.magmax,
                                         dismin=self.para.dismin, dismax=self.para.dismax)
            except Exception as e:
                raise ConnectionError(e)
        else:
            try:
                self.logger.RFlog.info(
                    'Searching earthquakes from {0} to {1}'.format(self.date_begin.strftime('%Y.%m.%dT%H:%M:%S'),
                                                                   self.date_end.strftime('%Y.%m.%dT%H:%M:%S')))
                self.eq_lst = read_catalog(self.para.catalogpath, self.para.date_begin, self.para.date_end,
                                           self.stainfo.stla, self.stainfo.stlo,
                                           magmin=self.para.magmin, magmax=self.para.magmax,
                                           dismin=self.para.dismin, dismax=self.para.dismax)
            except Exception as e:
                self.logger.RFlog.error('{0}'.format(e))
                raise e

    def match_eq(self, switchEN=False, reverseE=False, reverseN=False):
        try:
            self.logger.RFlog.info('Match SAC files')
            self.eqs = match_eq(self.eq_lst, self.para.datapath, self.stainfo.stla, self.stainfo.stlo,
                                ref_comp=self.para.ref_comp, suffix=self.para.suffix,
                                offset=self.para.offset, tolerance=self.para.tolerance,
                                dateformat=self.para.dateformat, switchEN=switchEN, reverseE=reverseE, reverseN=reverseN)
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise e
        self.logger.RFlog.info('{0} earthquakes matched'.format(self.eqs.shape[0]))

    def save(self, path=''):
        if path == '':
            path = '{0}.{1}.h5'.format(self.stainfo.network, self.stainfo.station)

        d = {'para': self.para.__dict__, 'stainfo': self.stainfo.__dict__, 'eq_lst': self.eq_lst, 'eqs': self.eqs}
        try:
            self.logger.RFlog.info('Saving project to {0}'.format(path))
            dd.io.save(path, d)
        except PerformanceWarning:
            pass
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise IOError(e)

    def load(self, path):
        try:
            self.logger.RFlog.info('Loading {0}'.format(path))
            fdd = dd.io.load(path)
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise IOError('Cannot read {0}'.format(path))

        try:
            self.para.__dict__.update(fdd['para'])
            self.stainfo.__dict__.update(fdd['stainfo'])
            self.eq_lst = fdd['eq_lst']
            self.eqs = fdd['eqs']
        except Exception as e:
            raise ValueError(e)

    def detrend(self):
        self.logger.RFlog.info('Detrend all data')
        for _, row in self.eqs.iterrows():
            row['data'].detrend()

    def filter(self, freqmin=0.05, freqmax=1., order=4):
        self.logger.RFlog.info('Filter all data from {0} to {1}'.format(freqmin, freqmax))
        for _, row in self.eqs.iterrows():
            row['data'].filter(freqmin=freqmin, freqmax=freqmax, order=order)

    def cal_phase(self):
        self.logger.RFlog.info('Calculate arrivals and ray parameters for all data')
        for _, row in self.eqs.iterrows():
            row['data'].get_arrival(self.model, row['evdp'], row['dis'])
            # row['data'].get_raypara(self.model, row['evdp'], row['dis'])

    def rotate(self, method='NE->RT'):
        self.logger.RFlog.info('Rotate {0} phase {1}'.format(self.para.phase, method))
        for _, row in self.eqs.iterrows():
            row['data'].rotate(row['bazi'], method=method, phase=self.para.phase)

    def drop_eq_snr(self, length=50):
        self.logger.RFlog.info('Reject data record with SNR less than {0}'.format(self.para.noisegate))
        drop_lst = []
        for i, row in self.eqs.iterrows():
            snr_R = row['data'].snr(length=length)
            if snr_R < self.para.noisegate:
                drop_lst.append(i)
        self.eqs.drop(drop_lst, inplace=True)
        self.logger.RFlog.info('{0} events left after SNR calculation'.format(self.eqs.shape[0]))

    def trim(self):
        for _, row in self.eqs.iterrows():
            row['data'].trim(self.para.time_before, self.para.time_after, self.para.phase)

    def deconv(self, criterion='crust', itmax=400, minderr=0.001):
        drop_lst = []
        if self.para.phase == 'P':
            shift = self.para.time_before
        elif self.para.phase == 'S':
            shift = self.para.time_after
            criterion = None
        else:
            pass

        for i, row in self.eqs.iterrows():
            row['data'].deconvolute(shift, self.para.time_after, self.para.gauss, phase=self.para.phase,
                                    only_r=self.para.only_r, itmax=itmax, minderr=minderr, target_dt=self.para.target_dt)
            if not row['data'].judge_rf(self.para.time_before, criterion=criterion):
                drop_lst.append(i)
                continue
            else:
                self.logger.RFlog.info('Iterative Decon {0} iterations: {1};'
                                       ' final RMS: {2}'.format(row['data'].datastr, row['data'].it,
                                                                row['data'].rms[-1]))
                row['data'].saverf(self.para.rfpath, phase=self.para.phase, shift=shift,
                                   evla=row['evla'], evlo=row['evlo'], evdp=row['evdp'], baz=row['bazi'], mag=row['mag'],
                                   gcarc=row['dis'], gauss=self.para.gauss, only_r=self.para.only_r)
        self.eqs.drop(drop_lst, inplace=True)


def InitRfProj(cfg_path):
    pjt = rf()
    pjt.para = CfgParser(cfg_path)
    return pjt


def rf_test():
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    logpath = '/Users/xumj/Codes/seispy/Scripts/EventCMT.dat'
    # logpath = '/home/xu_mijian/Codes/seispy/Scripts/EventCMT.dat'
    datapath = '/Users/xumj/Researches/test4seispy/data'
    # datapath = '/home/xu_mijian/xu_mijian/NJ2_SRF/data'
    proj_file = '/Users/xumj/Researches/test4seispy/test.h5'
    # proj_file = '/home/xu_mijian/xu_mijian/NJ2_SRF/test.h5'
    RFpath = '/Users/xumj/Researches/test4seispy/RFresult'
    # RFpath = '/home/xu_mijian/xu_mijian/NJ2_SRF/RFresult'

    rfproj = rf()
    # rfproj.load(proj_file)
    #
    # '''
    rfproj.date_begin = date_begin
    rfproj.date_end = date_end
    rfproj.para.datapath = datapath
    rfproj.para.catalogpath = logpath
    rfproj.para.rfpath = RFpath
    rfproj.para.gauss = 2.
    rfproj.para.offset = 0
    rfproj.load_stainfo()
    rfproj.search_eq(local=True)
    rfproj.match_eq()
    # rfproj.detrend()
    # rfproj.filter(freqmin=0.03, freqmax=0.5)
    # rfproj.cal_phase()
    # rfproj.drop_eq_snr(length=50)
    # rfproj.save(proj_file)
    # rfproj.trim()
    # rfproj.rotate()
    # rfproj.save(proj_file)
    # '''

    # rfproj.deconv()
    # rfproj.save(proj_file)


def srf_test():
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    logpath = '/Users/xumj/Codes/seispy/Scripts/EventCMT.dat'
    # logpath = '/home/xu_mijian/Codes/seispy/Scripts/EventCMT.dat'
    datapath = '/Users/xumj/Researches/test4seispy/data'
    # datapath = '/home/xu_mijian/xu_mijian/NJ2_SRF/data'
    proj_file = '/Users/xumj/Researches/test4seispy/test.h5'
    # proj_file = '/home/xu_mijian/xu_mijian/NJ2_SRF/test.h5'
    RFpath = '/Users/xumj/Researches/test4seispy/RFresult'
    # RFpath = '/home/xu_mijian/xu_mijian/NJ2_SRF/RFresult'

    rfproj = rf()
    # rfproj.load(proj_file)
    #
    # '''
    rfproj.date_begin = date_begin
    rfproj.date_end = date_end
    rfproj.para.datapath = datapath
    rfproj.para.catalogpath = logpath
    rfproj.para.rfpath = RFpath
    rfproj.para.time_before = 100
    rfproj.para.time_after = 30
    rfproj.para.gauss = 1
    rfproj.para.only_r = True
    rfproj.para.phase = 'S'
    rfproj.load_stainfo()
    rfproj.search_eq(local=True)
    rfproj.match_eq()
    rfproj.detrend()
    rfproj.filter(freqmin=0.03, freqmax=0.5)
    rfproj.cal_phase()
    rfproj.drop_eq_snr(length=50)
    rfproj.trim()
    rfproj.rotate(method='ZNE->LQT')
    rfproj.save(proj_file)
    # '''

    rfproj.deconv()


def get_events_test():
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    print(get_events(date_begin, date_end, 32.051701, 118.8544))


def proj_test():
    cfg_file = '/workspace/seispy_test/RF.cfg'
    pjt = rf(cfg_file)
    pjt.load_stainfo()
    pjt.search_eq(local=True)
    pjt.match_eq()
    pjt.cal_phase()
    for _, row in pjt.eqs.iterrows():
        row['data'].snr()


if __name__ == '__main__':
    # get_events_test()
    rf_test()
    # proj_test()
