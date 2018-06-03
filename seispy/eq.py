import numpy as np
import obspy
from obspy.io.sac import SACTrace
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


class eq(object):
    def __init__(self, pathname, datestr, suffix):
        self.datastr = datestr
        self.st = obspy.read(join(pathname, '*' + datestr + '.*.' + suffix))
        if not (self.st[0].stats.npts == self.st[1].stats.npts == self.st[2].stats.npts):
            raise ValueError('Samples are different in 3 components')
        self.st.sort()
        self.rf = obspy.Stream()
        self.timeoffset = 0
        self.rms = np.array([0])
        self.it = 0
        self.PArrival = None
        self.PRaypara = None
        self.SArrival = None
        self.SRaypara = None

    def __str__(self):
        return('Event data class {0}'.format(self.datastr))

    def detrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=order, zerophase=True)

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

        if self.rf == obspy.Stream():
            raise ValueError('Please cut out 3 components before')

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
            self.rf.rotate('->ZNE', inventory=inv)

        if method == 'NE->RT':
            self.rf.rotate('NE->RT', back_azimuth=bazi)
        elif method == 'RT->NE':
            self.rf.rotate('RT->NE', back_azimuth=bazi)
        elif method == 'ZNE->LQT':
            self.rf.rotate('ZNE->LQT', back_azimuth=bazi, inclination=inc)
        elif method == 'LQT->ZNE':
            self.rf.rotate('LQT->ZNE', back_azimuth=bazi, inclination=inc)
        else:
            pass

    def snr(self, length=50, phase='P'):
        st_noise = self.trim(length, 0, phase=phase, isreturn=True)
        st_signal = self.trim(0, length, phase=phase, isreturn=True)
        snr_R = seispy.geo.snr(st_signal[1].data, st_noise[1].data)
        return snr_R
    
    def get_time_offset(self, event_time):
        if not isinstance(event_time, obspy.core.utcdatetime.UTCDateTime):
            raise TypeError('Event time should be UTCDateTime type in obspy')
        self.timeoffset = self.st[2].stats.starttime - event_time

    def arr_correct(self, write_to_sac=True):
        """
        offset = sac.b - real o
        """
        Parr_time = self.PArrival.time
        Sarr_time = self.SArrival.time

        Pcorrect_time = Parr_time - self.timeoffset
        Scorrect_time = Sarr_time - self.timeoffset

        if write_to_sac:
            for tr in self.st:
                tr.stats.sac.t0 = Pcorrect_time
                tr.stats.sac.kt0 = 'P'
                tr.stats.sac.t1 = Scorrect_time
                tr.stats.sac.kt1 = 'S'

        return Pcorrect_time, Scorrect_time

    def trim(self, time_before, time_after, phase='P', isreturn=False):
        """
        offset = sac.b - real o
        """
        if phase not in ['P', 'S']:
            raise ValueError('Phase must in \'P\' or \'S\'')
        P_arr, S_arr = self.arr_correct(write_to_sac=False)
        time_dict = dict(zip(['P', 'S'], [P_arr, S_arr]))
        
        t1 = self.st[2].stats.starttime + (time_dict[phase] - time_before)
        t2 = self.st[2].stats.starttime + (time_dict[phase] + time_after)
        if isreturn:
            return self.st.copy().trim(t1, t2)
        else:
            self.rf = self.st.copy().trim(t1, t2)
    
    def deconvolute(self, shift, f0, phase='P', only_r=False, itmax=400, minderr=0.001):
        if phase not in ['P', 'S']:
            raise ValueError('Phase must in \'P\' or \'S\'')
        if self.rf == obspy.Stream():
            raise ValueError('Please run eq.trim at first')
        if phase == 'P':
            self.rf[1].data, self.rms, self.it = seispy.decov.decovit(self.rf[1].data, self.rf[2].data, self.rf[1].stats.delta,
                                 self.rf[1].data.shape[0], shift, f0, itmax, minderr)
            if not only_r:
                self.rf[0].data, self.rms, self.it = seispy.decov.decovit(self.rf[0].data, self.rf[2].data, self.rf[1].stats.delta,
                                     self.rf[1].data.shape[0], shift, f0, itmax, minderr)
        elif phase == 'S':
            if 'Q' not in self.rf[1].stats.channel or 'L' not in self.rf[2].stats.channel:
                raise ValueError('Please rotate component to \'LQT\'')
            Q = np.flip(self.rf[1].data, axis=0)
            L = np.flip(self.rf[2].data, axis=0)
            self.rf[2].data,self.rms, self.it = seispy.decov.decovit(L, -Q, self.rf[1].stats.delta, self.rf[1].shape[0], shift,
                                                         f0, itmax, minderr)
        else:
            pass

    def saverf(self, path, only_r=False):
        if only_r:
            loop_lst = (1)
        else:
            loop_lst = (0, 1)

        filename = join(path, self.datastr)
        for i in loop_lst:
            tr = SACTrace.from_obspy_trace(self.rf[i])
            tr.b = 0
            tr.write(filename + '_{0}.sac'.format(tr.kcmpnm[-1]))

    def judge_rf(self, shift, criterion='crust'):
        if not isinstance(criterion, str):
            raise TypeError('criterion should be string in [\'crust\', \'mtz\', \'lab\']')
        elif criterion.lower() not in ['crust', 'mtz', 'lab']:
            raise ValueError('criterion should be string in [\'crust\', \'mtz\', \'lab\']')
        else:
            criterion = criterion.lower()

        if criterion == 'crust':
            time_P1 = int(np.floor((-2+shift)/self.rf[1].stats.delta))
            time_P2 = int(np.floor((2+shift)/self.rf[1].stats.delta))
            max_P = np.max(self.rf[1].data[time_P1:time_P2])
            if max_P == np.max(np.abs(self.rf[1].data)) and max_P < 1:
                return True
            else:
                return False
        elif criterion == 'mtz':
            max_deep = np.max(np.abs(self.rf[1].data[int((25 + shift) / self.rf[1].stats.delta):]))
            time_P1 = int(np.floor((-2 + shift) / self.rf[1].stats.delta))
            time_P2 = int(np.floor((2 + shift) / self.rf[1].stats.delta))
            max_P = np.max(self.rf[1].data[time_P1:time_P2])
            if self.rms[-1] < 0.2 and max_deep < max_P * 0.3 and \
                  max_P == np.max(np.abs(self.rf[1].data)) and max_P < 1:
                return True
            else:
                return False
        elif criterion == 'lab':
            pass
        else:
            pass


if __name__ == '__main__':
    pass
