import obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.io.sac import SACTrace
from obspy.core.event.catalog import read_events
import re
from os.path import join, exists
from os import makedirs
from seispy.io import Query, _cat2df
from seispy.para import RFPara
from seispy import distaz
from seispy.eq import EQ
from seispy.setuplog import setuplog
from seispy.catalog import read_catalog_file
import glob
import numpy as np
from datetime import timedelta
import pandas as pd
import configparser
import argparse
import sys
import pickle


def pickphase(eqs, para, logger):
    from PySide6.QtWidgets import QApplication
    from seispy.pickseis.sviewerui import MatplotlibWidget
    app = QApplication(sys.argv)
    ui = MatplotlibWidget(eqs, para, logger)
    ui.show()
    if app.exec() == 0:
        ui.exit_app()
        return


def _add_header(st, evt_time, stainfo):
    sta = stainfo.query.stations[0][0]
    for tr in st:
        header = SACTrace.from_obspy_trace(tr).to_obspy_trace().stats.sac
        channel = tr.stats.channel
        for ch in sta.channels:
            if ch.code == channel:
                header.cmpaz = ch.azimuth
                header.cmpinc = -ch.dip
        header.stla = stainfo.stla
        header.stlo = stainfo.stlo
        header.stel = stainfo.stel
        header.b = 0
        header.o = evt_time - tr.stats.starttime
        tr.stats.__dict__['sac'] = header


class SACFileNotFoundError(Exception):
    def __init__(self, matchkey):
        self.matchkey = matchkey
    def __str__(self):
        print('No sac files found with {}'.format(self.matchkey))


def datestr2regex(datestr):
    pattern = datestr.replace('%y', r'\d{2}')
    pattern = pattern.replace('%Y', r'\d{4}')
    pattern = pattern.replace('%m', r'\d{2}')
    pattern = pattern.replace('%d', r'\d{2}')
    pattern = pattern.replace('%j', r'\d{3}')
    pattern = pattern.replace('%H', r'\d{2}')
    pattern = pattern.replace('%M', r'\d{2}')
    pattern = pattern.replace('%S', r'\d{2}')
    return pattern


def read_catalog(logpath:str, b_time, e_time, stla:float, stlo:float,
                 magmin=5.5, magmax=10., dismin=30., dismax=90.):
    """Read local catalog with seispy or QUAKEML format

    :param logpath: Path to catalogs
    :type logpath: str
    :param b_time: Start time 
    :type b_time: obspy.UTCDateTime
    :param e_time: End time
    :type e_time: obspy.UTCDateTime
    :param stla: Station latitude
    :type stla: float
    :param stlo: Station longitude
    :type stlo: float
    :param magmin: Minimum magnitude, defaults to 5.5
    :type magmin: float, optional
    :param magmax: Maximum magnitude, defaults to 10
    :type magmax: float, optional
    :param dismin: Minimum distance, defaults to 30
    :type dismin: float, optional
    :param dismax: Maximum distance, defaults to 90
    :type dismax: float, optional
    :return: list of earthquakes
    :rtype: pandas.DataFrame
    """
    try:
        eq_lst = read_catalog_file(logpath)
    except:
        events = read_events(logpath, 'QUAKEML')
        eq_lst = _cat2df(events)
    dis = distaz(stla, stlo, eq_lst['evla'], eq_lst['evlo']).delta
    eq_lst = eq_lst[(eq_lst['date']>=b_time) & (eq_lst['date']<=e_time) & \
                    (eq_lst['mag']>=magmin) & (eq_lst['mag']<=magmax) & \
                    (dis>=dismin) & (dis<=dismax)]
    return eq_lst


def fetch_waveform(eq_lst, para, model, logger):
    tb = np.max([2*para.noiselen, 2*para.time_before])
    te = np.max([2*para.noiselen, 2*para.time_after])
    query = Query(para.data_server)
    new_col = ['dis', 'bazi', 'data', 'datestr']
    eq_match = pd.DataFrame(columns=new_col)
    for i, row in eq_lst.iterrows():
        datestr = row['date'].strftime('%Y.%j.%H.%M.%S')
        daz = distaz(para.stainfo.stla, para.stainfo.stlo, row['evla'], row['evlo'])
        arrivals = model.get_travel_times(row['evdp'], daz.delta, phase_list=[para.phase])
        if not arrivals:
            logger.RFlog.error('The phase of {} with source depth {} and distance {} is not exists'.format(
                                     para.phase, row['evdp'], daz.delta))
            continue
        if len(arrivals) > 1:
            logger.RFlog.error('More than one phase were calculated with source depth of {} and distance of {}'.format(
                               row['evdp'], daz.delta))
        else:
            arr_time = arrivals[0].time
        t1 = row['date']+arr_time-tb
        t2 = row['date']+arr_time+te
        try:
            logger.RFlog.info('Fetch waveforms of ({}/{}) event {} from {}'.format(
                              i+1, eq_lst.shape[0], datestr, para.data_server))
            st = query.client.get_waveforms(para.stainfo.network, para.stainfo.station,
                                            para.stainfo.location, para.stainfo.channel, t1, t2)
            _add_header(st, row['date'], para.stainfo)
        except Exception as e:
            logger.RFlog.error('Error in fetching waveforms of event {}: {}'.format(datestr, str(e).strip()))
            continue
        try:
            this_eq = EQ.from_stream(st)
        except Exception as e:
            logger.RFlog.error('{}'.format(e))
            continue
        this_eq.get_time_offset(row['date'])
        this_df = pd.DataFrame([[daz.delta, daz.baz, this_eq, datestr]], columns=new_col, index=[i])
        eq_match = pd.concat([eq_match, this_df])
    ind = eq_match.index.drop_duplicates(keep=False)
    eq_match = eq_match.loc[ind]
    return pd.concat([eq_lst, eq_match], axis=1, join='inner')


def match_eq(eq_lst, pathname, stla, stlo, logger, ref_comp='Z', suffix='SAC', offset=None,
             tolerance=210, dateformat='%Y.%j.%H.%M.%S'):
    pattern = datestr2regex(dateformat)
    ref_eqs = glob.glob(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))
    if len(ref_eqs) == 0:
        raise SACFileNotFoundError(join(pathname, '*{0}*{1}'.format(ref_comp, suffix)))
    sac_files = []
    for ref_sac in ref_eqs:
        try:
            datestr = re.findall(pattern, ref_sac)[0]
        except IndexError:
            raise IndexError('Error data format of {} in {}'.format(pattern, ref_sac))
        if isinstance(offset, (int, float)):
            sac_files.append([datestr, UTCDateTime.strptime(datestr, dateformat), -offset])
        elif offset is None:
            try:
                tr = obspy.read(ref_sac)[0]
            except TypeError:
                continue
            sac_files.append([datestr, tr.stats.starttime-tr.stats.sac.b, tr.stats.sac.o])
        else:
            raise TypeError('offset should be int or float type')
    new_col = ['dis', 'bazi', 'data', 'datestr']
    eq_match = pd.DataFrame(columns=new_col)
    for datestr, b_time, offs in sac_files:
        date_range_begin = b_time + timedelta(seconds=offs - tolerance)
        date_range_end = b_time + timedelta(seconds=offs + tolerance)
        results = eq_lst[(eq_lst.date > date_range_begin) & (eq_lst.date < date_range_end)]
        if len(results) != 1:
            continue
        try:
            this_eq = EQ(pathname, datestr, suffix)
        except Exception as e:
            logger.RFlog.error(''.format(e))
            continue
        this_eq.get_time_offset(results.iloc[0]['date'])
        daz = distaz(stla, stlo, results.iloc[0]['evla'], results.iloc[0]['evlo'])
        this_df = pd.DataFrame([[daz.delta, daz.baz, this_eq, datestr]], columns=new_col, index=results.index.values)
        eq_match = pd.concat([eq_match, this_df])
    ind = eq_match.index.drop_duplicates(keep=False)
    eq_match = eq_match.loc[ind]
    return pd.concat([eq_lst, eq_match], axis=1, join='inner')


def CfgModify(cfg_file, session, key, value):
    cf = configparser.RawConfigParser()
    try:
        cf.read(cfg_file)
    except Exception:
        raise FileNotFoundError('Cannot open configure file %s' % cfg_file)
    cf.set(session, key, value)
    cf.write(open(cfg_file, 'w'))


class RF(object):
    def __init__(self, cfg_file=None, log=None):
        if log is None:
            self.logger = setuplog()
        else:
            self.logger = log
        if cfg_file is None:
            self.para = RFPara()
        elif isinstance(cfg_file, str):
            if not exists(cfg_file):
                self.logger.RFlog.error('No such file of {}.'.format(cfg_file))
                sys.exit(1)
            self.para = RFPara.read_para(cfg_file)
        else:
            raise TypeError('cfg should be in \'str\' format rather than \'{0}\''.format(type(cfg_file)))
        if not isinstance(self.para, RFPara):
            raise TypeError('Input value should be class seispy.rf.para')
        self.eq_lst = pd.DataFrame()
        self.eqs = pd.DataFrame()
        self.model = TauPyModel(self.para.velmod)
        self.stainfo = self.para.stainfo
        self.baz_shift = 0

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
            if self.para.use_remote_data:
                self.logger.RFlog.info('Load station info of {}.{} from {} web-service'.format(
                    self.para.stainfo.network, self.para.stainfo.station, self.para.data_server))
                self.para.stainfo.get_station_from_ws(self.para.data_server)
                try:
                    self.para._check_date_range()
                except Exception as e:
                    self.logger.RFlog.error('{}'.format(e))
                    sys.exit(1)
            else:
                self.logger.RFlog.info('Load station info from {0}'.format(self.para.datapath))
                self.para.stainfo.load_stainfo(self.para.datapath, self.para.ref_comp, self.para.suffix)
        except Exception as e:
            self.logger.RFlog.error('Error in loading station info: {0}'.format(e))
            sys.exit(1)

    def search_eq(self, local=False, catalog=None):
        if not local:
            try:
                self.logger.RFlog.info('Searching earthquakes from {}'.format(self.para.cata_server))
                if self.para.date_end > UTCDateTime():
                    self.para.date_end = UTCDateTime()
                query = Query(self.para.cata_server)
                query.get_events(starttime=self.para.date_begin, endtime=self.para.date_end,
                                 latitude=self.para.stainfo.stla, longitude=self.para.stainfo.stlo,
                                 minmagnitude=self.para.magmin, maxmagnitude=self.para.magmax,
                                 minradius=self.para.dismin, maxradius=self.para.dismax, catalog=catalog)
                self.eq_lst = query.events
            except Exception as e:
                raise ConnectionError(e)
        else:
            try:
                self.logger.RFlog.info(
                    'Searching earthquakes from {0} to {1}'.format(self.date_begin.strftime('%Y.%m.%dT%H:%M:%S'),
                                                                   self.date_end.strftime('%Y.%m.%dT%H:%M:%S')))
                self.eq_lst = read_catalog(self.para.catalogpath, self.para.date_begin, self.para.date_end,
                                           self.para.stainfo.stla, self.para.stainfo.stlo,
                                           magmin=self.para.magmin, magmax=self.para.magmax,
                                           dismin=self.para.dismin, dismax=self.para.dismax)
            except Exception as e:
                self.logger.RFlog.error('{0}'.format(e))
                raise e
        self.logger.RFlog.info('{} earthquakes are found'.format(self.eq_lst.shape[0]))

    def match_eq(self):
        try:
            if self.para.use_remote_data:
                self.logger.RFlog.info('Fetch seismic data from {}'.format(self.para.data_server))
                self.eqs = fetch_waveform(self.eq_lst, self.para, self.model, self.logger)
            else:
                self.logger.RFlog.info('Associating SAC files with earthquakes')
                self.eqs = match_eq(self.eq_lst, self.para.datapath, self.para.stainfo.stla, 
                                    self.para.stainfo.stlo, self.logger,
                                    ref_comp=self.para.ref_comp, suffix=self.para.suffix,
                                    offset=self.para.offset, tolerance=self.para.tolerance,
                                    dateformat=self.para.dateformat)
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise e
        if self.eqs.shape[0] == 0:
            self.logger.RFlog.warning('No earthquakes associated, please check configurations.'.format(self.eqs.shape[0]))
            sys.exit(1)
        else:
            self.logger.RFlog.info('{0} earthquakes are associated'.format(self.eqs.shape[0]))

    def save_raw_data(self):
        if not exists(self.para.datapath):
            makedirs(self.para.datapath)
        for i, row in self.eqs.iterrows():
            row['data'].write(self.para.datapath, row['date'])

    def savepjt(self):
        eqs = self.eqs.copy()
        for _, row in eqs.iterrows():
            row['data'].cleanstream()
        try:
            self.logger.RFlog.info('Saving project to {0}'.format(self.para.pjtpath))
            with open(self.para.pjtpath, 'wb') as f:
                pickle.dump({'para': self.para, 'eqs': eqs}, f, -1)
        except Exception as e:
            self.logger.RFlog.error('{0}'.format(e))
            raise IOError(e)
 
    @classmethod
    def loadpjt(cls, path):
        with open(path, 'rb') as f:
            rfdata = pickle.load(f)
        pjt = cls(phase=rfdata['para'].phase)
        pjt.para = rfdata['para']
        if not exists(pjt.para.datapath):
            pjt.logger.RFlog.error('Data path {} was not found'.format(pjt.para.datapath))
            sys.exit(1)
        pjt.load_stainfo()
        pjt.eqs = rfdata['eqs']
        for _, row in pjt.eqs.iterrows():
            try:
                row['data'].readstream()
            except Exception as e:
                pjt.logger.RFlog.warning('Cannot read {}, skipping'.format(row['data'].filestr))
                continue
        return pjt

    def channel_correct(self):
        if self.para.switchEN or self.para.reverseN or self.para.reverseE:
            self.logger.RFlog.info('Correct components with switchEN: {}, reverseE: {}, reverseN: {}'.format(
                                    self.para.switchEN, self.para.reverseE, self.para.reverseN))
            for _, row in self.eqs.iterrows():
                row['data'].channel_correct(self.para.switchEN, self.para.reverseE, self.para.reverseN)
        
    def detrend(self):
        self.logger.RFlog.info('Detrend all data')
        drop_idx = []
        for i, row in self.eqs.iterrows():
            try:
                row['data'].detrend()
            except:
                self.logger.RFlog.error('Data error in {}'.format(row['date'].strftime('%Y.%j.%H.%M.%S')))
                drop_idx.append(i)
        self.eqs.drop(drop_idx, inplace=True)

    def filter(self, freqmin=None, freqmax=None, order=4):
        if freqmin is None:
            freqmin = self.para.freqmin
        if freqmax is None:
            freqmax = self.para.freqmax
        self.logger.RFlog.info('Filter all data from {0} to {1}'.format(freqmin, freqmax))
        for _, row in self.eqs.iterrows():
            row['data'].filter(freqmin=freqmin, freqmax=freqmax, order=order)

    def cal_phase(self):
        self.logger.RFlog.info('Calculate {} arrivals and ray parameters for all data'.format(self.para.phase))
        for _, row in self.eqs.iterrows():
            row['data'].get_arrival(self.model, row['evdp'], row['dis'], phase=self.para.phase)

    def baz_correct(self, time_b=10, time_e=20, offset=90, correct_angle=None):
        if correct_angle is not None:
            self.logger.RFlog.info('correct back-azimuth with {} deg.'.format(correct_angle))
            self.eqs['bazi'] = np.mod(self.eqs['bazi'] + correct_angle, 360)
        else:
            self.logger.RFlog.info('correct back-azimuth with T energy minimization')
            y = np.arange(self.eqs.shape[0])
            shift_all = np.array([])
            x = np.arange(-offset, offset)
            ampt_all = np.empty([0, x.shape[0]])
            for i, row in self.eqs.iterrows():
                curr_baz, ampt = row['data'].search_baz(row['bazi'], time_b=time_b, time_e=time_e, offset=offset)
                shift_all = np.append(shift_all, curr_baz)
                ampt_all = np.vstack((ampt_all, ampt))
            if None in shift_all:
                self.logger.RFlog.error('Range of searching bazi is too small.')
                sys.exit(1)
            self.baz_shift = np.mean(shift_all[np.where(np.logical_not(np.isnan(shift_all)))])
            self.logger.RFlog.info('Average {:.1f} deg offset in back-azimuth'.format(self.baz_shift))

    def rotate(self, search_inc=False):
        targ_comp = ''.join(sorted(self.para.comp.upper()))
        if targ_comp == 'RTZ':
            method='NE->RT'
        elif targ_comp == 'LQT':
            method='ZNE->LQT'
        else:
            raise ValueError('comp must be in RTZ or LQT.')
        self.logger.RFlog.info('Rotate {0} phase to {1}'.format(self.para.phase, method))
        drop_idx = []
        for i, row in self.eqs.iterrows():
            try:
                row['data'].rotate(row['bazi'], method=method, search_inc=search_inc, baz_shift=self.baz_shift)
            except Exception as e:
                self.logger.RFlog.error('{}: {}'.format(row['data'].datestr, e))
                drop_idx.append(i)
                continue
            if search_inc:
                self.logger.RFlog.info('The incidence angle of {} was corrected by {:.1f} deg'.format(
                                       row['data'].datestr, row['data'].inc_correction))
        self.eqs.drop(drop_idx, inplace=True)

    def drop_eq_snr(self, length=None):
        if length is None:
            length = self.para.noiselen
        self.logger.RFlog.info('Reject data record with SNR less than {0}'.format(self.para.noisegate))
        drop_lst = []
        for i, row in self.eqs.iterrows():
            snr_E, snr_N, snr_Z = row['data'].snr(length=length)
            mean_snr = np.mean([snr_E, snr_N, snr_Z])
            if mean_snr < self.para.noisegate:
                drop_lst.append(i)
        self.eqs.drop(drop_lst, inplace=True)
        self.logger.RFlog.info('{0} events left after SNR calculation'.format(self.eqs.shape[0]))

    def trim(self):
        self.logger.RFlog.info('Trim waveforms from {0:.2f} before {2} to {1:.2f} after {2}'.format(
                               self.para.time_before, self.para.time_after, self.para.phase))
        for _, row in self.eqs.iterrows():
            row['data'].trim(self.para.time_before, self.para.time_after)
    
    def pick(self, prepick=True, stl=5, ltl=10):
        if prepick:
            self.logger.RFlog.info('Pre-pick {} arrival using STA/LTA method'.format(self.para.phase))
            for _, row in self.eqs.iterrows():
                row['data'].phase_trigger(self.para.time_before, self.para.time_after,
                                          stl=stl, ltl=ltl)
        self.logger.RFlog.info('{0} events left after virtual checking'.format(self.eqs.shape[0]))
        pickphase(self.eqs, self.para, self.logger)

    def deconv(self):
        shift = self.para.time_before
        time_after = self.para.time_after
        drop_lst = []

        count = 0
        for i, row in self.eqs.iterrows():
            count += 1
            try:
                row['data'].deconvolute(shift, time_after, method=self.para.decon_method, f0=self.para.gauss,
                                        only_r=self.para.only_r, itmax=self.para.itmax, minderr=self.para.minderr,
                                        wlevel=self.para.wlevel, target_dt=self.para.target_dt)
                if self.para.decon_method == 'iter':
                    self.logger.RFlog.info('Iterative Decon {0} ({3}/{4}) iterations: {1}; final RMS: {2:.4f}'.format(
                        row['data'].datestr, row['data'].rf[0].stats.iter,
                        row['data'].rf[0].stats.rms[-1], count, self.eqs.shape[0]))
                elif self.para.decon_method == 'water':
                    self.logger.RFlog.info('Water level Decon {} ({}/{}); RMS: {:.4f}'.format(
                        row['data'].datestr, count, self.eqs.shape[0], row['data'].rf[0].stats.rms))
            except Exception as e:
                self.logger.RFlog.error('{}: {}'.format(row['data'].datestr, e))
                drop_lst.append(i)
        self.eqs.drop(drop_lst, inplace=True)

    def saverf(self, gauss=None):
        npts = int((self.para.time_before + self.para.time_after)/self.para.target_dt+1)
        if self.para.phase[-1] == 'P':
            shift = self.para.time_before
        elif self.para.phase[-1] == 'S':
            shift = self.para.time_after
        else:
            pass
        good_lst = []
        if isinstance(self.para.gauss, (int, float)):
            gauss = self.para.gauss
        elif gauss is None:
            gauss = self.para.gauss[0]
        else:
            if gauss in self.para.gauss:
                pass
            else:
                raise ValueError('gauss should be a element in the seispy.para.RFPara.gauss')

        if self.para.rmsgate is not None:
            self.logger.RFlog.info('Save RFs with final RMS less than {:.2f} and criterion of {}'.format(self.para.rmsgate, self.para.criterion))
        else:
            self.logger.RFlog.info('Save RFs with and criterion of {}'.format(self.para.criterion))
        for i, row in self.eqs.iterrows():
            if row['data'].judge_rf(gauss, shift, npts, criterion=self.para.criterion, rmsgate=self.para.rmsgate):
                row['data'].saverf(self.para.rfpath, evtstr=row['date'].strftime('%Y.%j.%H.%M.%S'), shift=shift,
                                   evla=row['evla'], evlo=row['evlo'], evdp=row['evdp'], baz=row['bazi'],
                                   mag=row['mag'], gcarc=row['dis'], gauss=gauss, only_r=self.para.only_r,
                                   user9=self.baz_shift)
                good_lst.append(i)
        self.logger.RFlog.info('{} PRFs are saved.'.format(len(good_lst)))
        self.eqs = self.eqs.loc[good_lst]


def setpar():
    parser = argparse.ArgumentParser(description="Set parameters to configure file")
    parser.add_argument('cfg_file', type=str, help='Path to configure file')
    parser.add_argument('session', type=str, help='session name')
    parser.add_argument('key', type=str, help='key name')
    parser.add_argument('value', type=str, help='value')
    arg = parser.parse_args()
    CfgModify(arg.cfg_file, arg.session, arg.key, arg.value)


if __name__ == '__main__':
    pass
    # get_events_test()
    # prf()
    # proj_test()
