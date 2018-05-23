#!/usr/bin/env python

import numpy as np
import obspy
import re
import io
from os.path import dirname, join, expanduser
import seispy
from seispy.setuplog import setuplog
import glob
from datetime import timedelta, datetime
import pandas as pd
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
    ex_sac = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))[0]
    ex_tr = obspy.read(ex_sac)[0]
    return ex_tr.stats.network, ex_tr.stats.station, ex_tr.stats.sac.stla, ex_tr.stats.sac.stlo


def match_eq(eq_lst, pathname, stla, stlo, ref_comp='Z', suffix='SAC', offset=0,
             tolerance=210, dateformat='%Y.%j.%H.%M.%S'):
    pattern = datestr2regex(dateformat)
    ref_eqs = glob.glob(join(pathname, '*{0}*.{1}'.format(ref_comp, suffix)))
    sac_files = []
    for ref_sac in ref_eqs:
        datestr = re.findall(pattern, ref_sac)[0]
        tr = obspy.read(ref_sac)[0]
        sac_files.append([datestr, tr])
    new_col = ['dis', 'bazi', 'data']
    eq_match = pd.DataFrame(columns=new_col)
    for i, evt in eq_lst.iterrows():
        tmp_datestr = []
        for datestr, tr in sac_files:
            if tr.stats.starttime - timedelta(seconds=offset + tolerance) <= evt['date'] \
                    <= tr.stats.starttime + timedelta(seconds=-offset + tolerance):
                tmp_datestr.append(datestr)
        if len(tmp_datestr) == 1:
            this_eq = eq(pathname, tmp_datestr[0], suffix)
            daz = seispy.distaz(stla, stlo, evt['evla'], evt['evlo'])
            this_df = pd.DataFrame([[daz.delta, daz.baz, this_eq]], columns=new_col, index=[i])
            eq_match = eq_match.append(this_df)
    return pd.concat([eq_lst, eq_match], axis=1, join='inner')


class eq(object):
    def __init__(self, pathname, datestr, suffix):
        self._datastr = datestr
        self.st = obspy.read(join(pathname, '*' + datestr + '.*.' + suffix))
        self.rf = obspy.Stream()
        self.PArrival = None
        self.PRaypara = None
        self.SArrival = None
        self.SRaypara = None

    def __str__(self):
        return('Event data class {0}'.format(self._datastr))

    def detrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=order)

    def get_arrival(self, model, evdp, dis):
        arrivals = model.get_travel_times(evdp, dis, phase_list=['P', 'S'])
        self.PArrival = arrivals[0]
        self.SArrival = arrivals[1]

    def get_raypara(self, model, evdp, dis):
        raypara = model.get_ray_paths(evdp, dis, phase_list=['P', 'S'])
        self.PRaypara = raypara[0]
        self.SRaypara = raypara[1]

    def rotate(self, bazi, inc=None, method='NE->RT', phase='P', inv=None):
        if phase not in ('P', 'S'):
            raise ValueError('phase must be in [\'P\', \'S\']')

        if inc is None:
            if self.PRaypara is None or self.SRaypara is None:
                raise ValueError('inc must be specified')
            elif phase == 'P':
                inc = self.PRaypara.incident_angle
            elif phase == 'S':
                inc = self.SRaypara.incident_angle
            else:
                pass

        if inv is not None:
            self.st.rotate('->ZNE', inventory=inv)

        if method == 'NE->RT':
            self.st.rotate('NE->RT', back_azimuth=bazi)
        elif method == 'RT->NE':
            self.st.rotate('RT->NE', back_azimuth=bazi)
        elif method == 'ZNE->LQT':
            self.st.rotate('ZNE->LQT', back_azimuth=bazi, inclination=inc)
        elif method == 'LQT->ZNE':
            self.st.rotate('LQT->ZNE', back_azimuth=bazi, inclination=inc)
        else:
            pass

    def snr(self, length=50):
        pass
    
    def time_correct(self, offset=0, write_to_sac=True):
        """
        offset = sac.b - real o
        """
        Parr_time = self.PArrival.time
        Sarr_time = self.SArrival.time

        Pcorrect_time = Parr_time - offset
        Scorrect_time = Sarr_time - offset

        if write_to_sac:
            for tr in self.st:
                tr.stats.sac.t0 = Pcorrect_time
                tr.stats.sac.kt0 = 'P'
                tr.stats.sac.t1 = Scorrect_time
                tr.stats.sac.kt1 = 'S'

        return Pcorrect_time, Scorrect_time

    def trim(self, time_before, time_after, phase='P', offset=0):
        """
        offset = sac.b - real o
        """
        if phase not in ['P', 'S']:
            raise ValueError('Phase must in \'P\' or \'S\'')
        dt = self.st[0].stats.delta
        P_arr, S_arr = self.time_correct(offset=offset, write_to_sac=False)
        time_dict = dict(zip(['P', 'S'], [P_arr, S_arr]))
        
        t1 = self.st[2].stats.starttime + (time_dict[phase] - time_before)
        t2 = self.st[2].stats.starttime + (time_dict[phase] + time_after)
        return self.st.copy().trim(t1, t2)


class para():
    def __init__(self):
        self.datapath = expanduser('~')
        self.RFpath = expanduser('~')
        self.imagepath = expanduser('~')
        self.catalogpath = join(expanduser('~'), 'EventCMT.dat')
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
        self.noiseget = 7
        self.gauss = 2
        self.target_dt = 0.01

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


def CfgParser(cfg_file):
    cf = configparser.ConfigParser()
    pa = para()
    try:
        cf.read(cfg_file)
    except Exception:
        raise FileNotFoundError('Cannot open configure file %s' % cfg_file)

    for key, value in cf.items('path'):
        pa.__dict__[key] = value
    for key, value in cf.items('para'):
        if key == 'date_begin':
            pa.__dict__[key] = obspy.UTCDateTime(value)
        elif key == 'date_end':
            pa.__dict__[key] = obspy.UTCDateTime(value)
        else:
            try:
                pa.__dict__[key] = float(value)
            except ValueError:
                pa.__dict__[key] = value
    return pa


class rf(object):
    def __init__(self):
        self.para = para()
        self.stainfo = stainfo()
        self.eq_lst = pd.DataFrame()
        self.eqs = pd.DataFrame()
        self.model = TauPyModel('iasp91')
        self.logger = setuplog()

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

    def match_eq(self):
        try:
            self.logger.RFlog.info('Match SAC files')
            self.eqs = match_eq(self.eq_lst, self.para.datapath, self.stainfo.stla, self.stainfo.stlo,
                                ref_comp=self.para.ref_comp, suffix=self.para.suffix,
                                offset=self.para.offset, tolerance=self.para.tolerance,
                                dateformat=self.para.dateformat)
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

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.logger.RFlog.info('Filter all data from {0} to {1}'.format(freqmin, freqmax))
        for _, row in self.eqs.iterrows():
            row['data'].filter(freqmin=freqmin, freqmax=freqmax, order=order)

    def phase(self):
        self.logger.RFlog.info('Calculate arrivals and ray parameters for all data')
        for _, row in self.eqs.iterrows():
            row['data'].get_arrival(self.model, row['evdp'], row['dis'])
            row['data'].get_raypara(self.model, row['evdp'], row['dis'])

    def rotate(self, method='NE->RT', phase='P'):
        self.logger.RFlog.info('Rotate {0} phase {1}'.format(phase, method))
        for _, row in self.eqs.iterrows():
            row['data'].rotate(row['bazi'], method='NE->RT', phase='P')


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
    # proj_file = '/home/xu_mijian/xu_mijian/NJ2_SRF/test.h5'
    proj_file = '/Users/xumj/Researches/test4seispy/test.h5'

    rfproj = rf()
    # rfproj.load(proj_file)
    rfproj.date_begin = date_begin
    rfproj.date_end = date_end
    rfproj.datapath = datapath
    rfproj.para.catalogpath = logpath
    rfproj.load_stainfo()
    rfproj.search_eq()
    rfproj.match_eq()
    # rfproj.save(proj_file)
    # rfproj.detrend()
    # rfproj.filter()
    # rfproj.phase()
    # rfproj.rotate()

    print(rfproj.eqs)


def get_events_test():
    date_begin = obspy.UTCDateTime('20130101')
    date_end = obspy.UTCDateTime('20140101')
    print(get_events(date_begin, date_end, 32.051701, 118.8544))


if __name__ == '__main__':
    # get_events_test()
    rf_test()
