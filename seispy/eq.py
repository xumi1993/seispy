import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.rotate import rotate2zne, rotate_zne_lqt
from scipy.signal import resample
from os.path import dirname, join, expanduser
import seispy


def rotateZNE(st):
    try:
        zne = rotate2zne(
            st[0], st[0].stats.sac.cmpaz, st[0].stats.sac.cmpinc,
            st[1], st[1].stats.sac.cmpaz, st[1].stats.sac.cmpinc,
            st[2], st[2].stats.sac.cmpaz, st[2].stats.sac.cmpinc)
    except Exception as e:
        raise ValueError('No specified cmpaz or cmpinc. {}'.format(e))
    for tr, new_data, component in zip(st, zne, "ZNE"):
        tr.data = new_data
        tr.stats.channel = tr.stats.channel[:-1] + component


class eq(object):
    def __init__(self, pathname, datestr, suffix, switchEN=False, reverseE=False, reverseN=False):
        self.datastr = datestr
        self.st = obspy.read(join(pathname, '*' + datestr + '*' + suffix))
        if not (self.st[0].stats.npts == self.st[1].stats.npts == self.st[2].stats.npts):
            raise ValueError('Samples are different in 3 components')
        if reverseE:
            self.st.select(channel='*[E2]')[0].data *= -1
        if reverseN:
            self.st.select(channel='*[N1]')[0].data *= -1
        if switchEN:
            chE = self.st.select(channel='*[E2]')[0].stats.channel
            chN = self.st.select(channel='*[N1]')[0].stats.channel
            self.st.select(channel='*[E2]')[0].stats.channel = chN
            self.st.select(channel='*[N1]')[0].stats.channel = chE
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

    def search_inc(self, bazi):
        inc_range = np.arange(0.1, 90, 0.1)
        # bazi_range = np.repeat(bazi, len(inc_range))
        # M_all = seispy.geo.rot3D(bazi=bazi_range, inc=inc_range)
        # ZEN = np.array([self.rf[2].data, self.rf[0].data, self.rf[1].data])
        s_range = self.trim(10, 10, phase='S', isreturn=True)
        # LQT_all = np.zeros([ZEN.shape[0], ZEN.shape[1], M_all.shape[2]])
        power = np.zeros(inc_range.shape[0])
        for i in range(len(inc_range)):
            # LQT_all[:, :, i] = M_all[:, :, i].dot(ZEN)
            l_comp, _, _ = rotate_zne_lqt(s_range[2].data, s_range[1].data, s_range[0].data, bazi, inc_range[i])
            power[i] = np.mean(l_comp ** 2)

        # real_inc_idx = seispy.geo.extrema(power, opt='min')
        real_inc_idx = np.argmin(power)
        # print(real_inc_idx)
        real_inc = inc_range[real_inc_idx]
        return real_inc

    def rotate(self, bazi, inc=None, method='NE->RT', phase='P', inv=None):
        if phase not in ('P', 'S'):
            raise ValueError('phase must be in [\'P\', \'S\']')

        if self.rf == obspy.Stream():
            raise ValueError('Please cut out 3 components before')

        if inc is None:
            if self.PArrival is None or self.SArrival is None:
                raise ValueError('inc must be specified')
            elif phase == 'P':
                inc = self.PArrival.incident_angle
            elif phase == 'S':
                inc = self.SArrival.incident_angle
                # inc = self.search_inc(bazi)
            else:
                pass

        if inv is not None:
            self.rf.rotate('->ZNE', inventory=inv)
        elif self.rf.select(channel='*1'):
            if self.rf.select(channel='*1')[0].stats.sac.cmpaz == 0:
                self.rf.select(channel='*1')[0].stats.channel = self.rf.select(channel='*1')[0].stats.channel[:-1] + 'N'
                self.rf.select(channel='*2')[0].stats.channel = self.rf.select(channel='*2')[0].stats.channel[:-1] + 'E'
            elif self.rf.select(channel='*1')[0].stats.sac.cmpaz != 0:
                rotateZNE(self.rf)
            self.rf.sort(['channel'])
        else:
            pass

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
        try:
            snr_R = seispy.geo.snr(st_signal[1].data, st_noise[1].data)
        except IndexError:
            snr_R = 0
        return snr_R
    
    def get_time_offset(self, event_time):
        if not isinstance(event_time, obspy.core.utcdatetime.UTCDateTime):
            raise TypeError('Event time should be UTCDateTime type in obspy')
        self.timeoffset = self.st[2].stats.starttime - event_time

    def arr_correct(self, write_to_sac=True):
        """
        offset = sac.b - real o
        """
        try:
            Parr_time = self.PArrival.time
            Sarr_time = self.SArrival.time
        except AttributeError:
            raise Exception('Please calculate PArrival or SArrival first')

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
    
    def deconvolute(self, shift, time_after, f0, phase='P', only_r=False, itmax=400, minderr=0.001, target_dt=None):
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
        else:
            if 'Q' not in self.rf[1].stats.channel or 'L' not in self.rf[2].stats.channel:
                raise ValueError('Please rotate component to \'LQT\'')
            Q = np.flip(self.rf[1].data, axis=0)
            L = np.flip(self.rf[2].data, axis=0)
            self.rf[2].data, self.rms, self.it = seispy.decov.decovit(L, -Q, self.rf[1].stats.delta, self.rf[1].data.shape[0], shift,
                                                         f0, itmax, minderr)
        if target_dt is not None:
            if self.rf[0].stats.delta != target_dt:
                # self.rf.resample(1 / target_dt)
                for tr in self.rf:
                    # tr.data = tr.data[0:-1]
                    tr.data = resample(tr.data, int((shift + time_after)/target_dt+1))

    def saverf(self, path, phase='P', shift=0, time_after=120, evla=-12345., evlo=-12345., evdp=-12345., mag=-12345.,
               gauss=0, baz=-12345., gcarc=-12345., only_r=False, **kwargs):
        if phase == 'P':
            if only_r:
                loop_lst = [1]
            else:
                loop_lst = [0, 1]
            rayp = seispy.geo.srad2skm(self.PArrival.ray_param)
        elif phase == 'S':
            loop_lst = [2]
            rayp = seispy.geo.srad2skm(self.SArrival.ray_param)
        else:
            raise ValueError('Phase must be in \'P\' or \'S\'')

        filename = join(path, self.datastr)
        for i in loop_lst:
            header = {'evla': evla, 'evlo': evlo, 'evdp': evdp, 'mag': mag, 'baz': baz,
                      'gcarc': gcarc, 'user0': rayp, 'kuser0': 'Ray Para', 'user1': gauss, 'kuser1': 'G factor'}
            for key in kwargs:
                header[key] = kwargs[key]
            for key, value in header.items():
                self.rf[i].stats['sac'][key] = value
            tr = SACTrace.from_obspy_trace(self.rf[i])
            tr.b = -shift
            tr.e = time_after
            tr.o = 0
            tr.write(filename + '_{0}_{1}.sac'.format(phase, tr.kcmpnm[-1]))

    def judge_rf(self, shift, criterion='crust'):
        if not isinstance(criterion, (str, type(None))):
            raise TypeError('criterion should be string in [\'crust\', \'mtz\']')
        elif criterion is None:
            pass
        elif criterion.lower() not in ['crust', 'mtz', 'lab']:
            raise ValueError('criterion should be string in [\'crust\', \'mtz\']')
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
            max_deep = np.max(np.abs(self.rf[1].data[int((30 + shift) / self.rf[1].stats.delta):]))
            time_P1 = int(np.floor((-5 + shift) / self.rf[1].stats.delta))
            time_P2 = int(np.floor((5 + shift) / self.rf[1].stats.delta))
            max_P = np.max(self.rf[1].data[time_P1:time_P2])
            if self.rms[-1] < 0.4 and max_deep < max_P * 0.4 and \
                  max_P == np.max(np.abs(self.rf[1].data)) and max_P < 1:
                return True
            else:
                return False
        elif criterion is None:
            return True
        else:
            pass


if __name__ == '__main__':
    pass
