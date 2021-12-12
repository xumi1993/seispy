import numpy as np
from numpy.lib.arraysetops import isin
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.rotate import rotate2zne, rotate_zne_lqt
from scipy.signal import resample
from os.path import dirname, join, expanduser
from seispy.decon import RFTrace
from seispy.geo import snr, srad2skm, rotateSeisENtoTR, skm2srad, \
                       rssq, extrema
from obspy.signal.trigger import recursive_sta_lta


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
    def __init__(self, pathname, datestr, suffix):
        self.datestr = datestr
        self.st = obspy.read(join(pathname, '*' + datestr + '*' + suffix))
        if len(self.st) < 3:
            channel = ' '.join([tr.stats.channel for tr in self.st])
            raise ValueError('Sismogram must be in 3 components, but there are only channel {} of {}'.format(channel, datestr))
        elif len(self.st) > 3:
            raise ValueError('{} has more than 3 components, please select to delete redundant seismic components'.format(datestr))
        else:
            pass
        # if not (self.st[0].stats.npts == self.st[1].stats.npts == self.st[2].stats.npts):
        #     raise ValueError('Samples are different in 3 components')
        self.st.sort()
        self.rf = obspy.Stream()
        self.timeoffset = 0
        self.rms = np.array([0])
        self.it = 0
        self.PArrival = None
        self.PRaypara = None
        self.SArrival = None
        self.SRaypara = None
        self.trigger_shift = 0

    def channel_correct(self, switchEN=False, reverseE=False, reverseN=False):
        if reverseE:
            self.st.select(channel='*[E2]')[0].data *= -1
        if reverseN:
            self.st.select(channel='*[N1]')[0].data *= -1
        if switchEN:
            chE = self.st.select(channel='*[E2]')[0].stats.channel
            chN = self.st.select(channel='*[N1]')[0].stats.channel
            self.st.select(channel='*[E2]')[0].stats.channel = chN
            self.st.select(channel='*[N1]')[0].stats.channel = chE

    def __str__(self):
        return('Event data class {0}'.format(self.datestr))

    def detrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')
        self.fix_channel_name()

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
        s_range = self.trim(20, 20, phase='S', isreturn=True)
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

    def search_baz(self, bazi, time_b=10, time_e=20, offset=90):
        p_arr, _ = self.arr_correct(write_to_sac=False)
        this_st = self.st.copy()
        this_st.filter('bandpass', freqmin=0.03, freqmax=0.5)
        this_st.trim(this_st[0].stats.starttime+p_arr-time_b, this_st[0].stats.starttime+p_arr+time_e)
        bazs = np.arange(-offset, offset) + bazi
        ampt = np.zeros_like(bazs)
        for i, b in enumerate(bazs):
            t, _ = rotateSeisENtoTR(this_st[0].data, this_st[1].data, b)
            ampt[i] = rssq(t)
        ampt = ampt / np.max(ampt)
        idx = extrema(ampt, opt='min')
        if len(idx) == 0:
            corr_baz = np.nan
        elif len(idx) > 1:
            corr_baz = None
        else:
            corr_baz = bazs[idx[0]] - bazi
        return corr_baz, ampt

    def fix_channel_name(self):
        if self.st.select(channel='??1') and self.st.select(channel='??Z') and hasattr(self.st.select(channel='*1')[0].stats.sac, 'cmpaz'):
            if self.st.select(channel='*1')[0].stats.sac.cmpaz == 0:
                self.st.select(channel='*1')[0].stats.channel = self.st.select(channel='*1')[0].stats.channel[:-1] + 'N'
                self.st.select(channel='*2')[0].stats.channel = self.st.select(channel='*2')[0].stats.channel[:-1] + 'E'
            elif self.st.select(channel='*1')[0].stats.sac.cmpaz != 0:
                rotateZNE(self.st)
            self.st.sort(['channel'])
        elif self.st.select(channel='*1'):
            self.st.select(channel='*1')[0].stats.channel = 'Z'
            self.st.select(channel='*2')[0].stats.channel = 'N'
            self.st.select(channel='*3')[0].stats.channel = 'E'
            self.st.sort(['channel'])
        else:
            pass

    def rotate(self, bazi, inc=None, method='NE->RT', phase='P', search_inc=False):
        if phase not in ('P', 'S'):
            raise ValueError('phase must be in [\'P\', \'S\']')

        if inc is None:
            if self.PArrival is None or self.SArrival is None:
                raise ValueError('inc must be specified')
            elif phase == 'P':
                inc = self.PArrival.incident_angle
            elif phase == 'S':
                if search_inc:
                    inc = self.search_inc(bazi)
                    #  if inc == 0.1 or inc == 89.9:
                        #  print(self.datestr)
                else:
                    inc = self.SArrival.incident_angle
                # inc = self.search_inc(bazi)
            else:
                pass

        if method == 'NE->RT':
            self.comp = 'rtz'
            self.st.rotate('NE->RT', back_azimuth=bazi)
        elif method == 'RT->NE':
            self.st.rotate('RT->NE', back_azimuth=bazi)
        elif method == 'ZNE->LQT':
            self.comp = 'lqt'
            self.st.rotate('ZNE->LQT', back_azimuth=bazi, inclination=inc)
        elif method == 'LQT->ZNE':
            self.st.rotate('LQT->ZNE', back_azimuth=bazi, inclination=inc)
        else:
            pass

    def snr(self, length=50, phase='P'):
        st_noise = self.trim(length, 0, phase=phase, isreturn=True)
        st_signal = self.trim(0, length, phase=phase, isreturn=True)
        try:
            snr_E = snr(st_signal[0].data, st_noise[0].data)
        except IndexError:
            snr_E = 0
        try:
            snr_N = snr(st_signal[1].data, st_noise[1].data)
        except IndexError:
            snr_N = 0
        try:
            snr_Z = snr(st_signal[2].data, st_noise[2].data)
        except IndexError:
            snr_Z = 0
        return snr_E, snr_N, snr_Z
    
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

    def _get_time(self, time_before, time_after, phase='P'):
        if phase not in ['P', 'S']:
            raise ValueError('Phase must in \'P\' or \'S\'')
        P_arr, S_arr = self.arr_correct(write_to_sac=False)
        time_dict = dict(zip(['P', 'S'], [P_arr, S_arr]))

        t1 = self.st[2].stats.starttime + (time_dict[phase] + self.trigger_shift - time_before)
        t2 = self.st[2].stats.starttime + (time_dict[phase] + self.trigger_shift + time_after)
        return t1, t2

    def phase_trigger(self, time_before, time_after, phase='S', stl=5, ltl=10):
        t1, t2 = self._get_time(time_before, time_after, phase)
        self.t1_pick, self.t2_pick = t1, t2
        self.st_pick = self.st.copy().trim(t1, t2)
        if phase == 'P':
            tr = self.st_pick.select(channel='*Z')[0]
        else:
            tr = self.st_pick.select(channel='*T')[0]
        df = tr.stats.sampling_rate
        cft = recursive_sta_lta(tr.data, int(stl*df), int(ltl*df))
        n_trigger = np.argmax(np.diff(cft)[int(ltl*df):])+int(ltl*df)
        self.t_trigger = t1 + n_trigger/df
        self.trigger_shift = n_trigger/df - time_before

    def trim(self, time_before, time_after, phase='P', isreturn=False):
        """
        offset = sac.b - real o
        """
        t1, t2 = self._get_time(time_before, time_after, phase)
        if isreturn:
            return self.st.copy().trim(t1, t2)
        else:
            self.st.trim(t1, t2)

    def deconvolute(self, shift, time_after, f0=2, phase='P', method='iter', only_r=False, itmax=400, minderr=0.001, wlevel=0.05, target_dt=None):
        self.method = method
        if phase not in ['P', 'S']:
            raise ValueError('Phase must in \'P\' or \'S\'')
        if method == 'iter':
            kwargs = {'method': method,
                      'f0': f0,
                      'tshift': shift,
                      'itmax': itmax,
                      'minderr': minderr}
        elif method == 'water':
            kwargs = {'method': method,
                      'f0': f0,
                      'tshift': shift,
                      'wlevel': wlevel}
        else:
            raise ValueError('method must be in \'iter\' or \'water\'')

        if phase == 'P':
            self.decon_p(**kwargs)
            if not only_r:
                self.decon_p(tcomp=True, **kwargs)
        else:
            # TODO: if 'Q' not in self.rf[1].stats.channel or 'L' not in self.rf[2].stats.channel:
            #     raise ValueError('Please rotate component to \'LQT\'')
            self.decon_s(**kwargs)
        if target_dt is not None:
            if self.rf[0].stats.delta != target_dt:
                # self.rf.resample(1 / target_dt)
                for tr in self.rf:
                    # tr.data = tr.data[0:-1]
                    tr.data = resample(tr.data, int((shift + time_after)/target_dt+1))
                    tr.stats.delta = target_dt
    
    def decon_p(self, tshift, tcomp=False, **kwargs):
        if self.comp == 'lqt':
            win = self.st.select(channel='*L')[0]
            if tcomp:
                uin = self.st.select(channel='*T')[0]
            else:
                uin = self.st.select(channel='*Q')[0]
        else:
            win = self.st.select(channel='*Z')[0]
            if tcomp:
                uin = self.st.select(channel='*T')[0]
            else:
                uin = self.st.select(channel='*R')[0]
        uout = RFTrace.deconvolute(uin, win, phase='P', tshift=tshift, **kwargs)
        self.rf.append(uout)


    def decon_s(self, tshift, **kwargs):
        if self.comp == 'lqt':
            win = self.st.select(channel='*Q')[0]
            uin = self.st.select(channel='*L')[0]
        else:
            win = self.st.select(channel='*R')[0]
            uin = self.st.select(channel='*Z')[0]
            win.data *= -1
        win.data[0:int((tshift-4)/win.stats.delta)] = 0
        uout = RFTrace.deconvolute(uin, win, phase='S', tshift=tshift, **kwargs)
        uout.data = np.flip(uout.data)
        self.rf.append(uout)


    def saverf(self, path, evtstr=None, phase='P', shift=0, evla=-12345., evlo=-12345., evdp=-12345., mag=-12345.,
               gauss=0, baz=-12345., gcarc=-12345., only_r=False, **kwargs):
        if phase == 'P':
            if self.comp == 'lqt':
                svcomp = 'Q'
            else:
                svcomp = 'R'
            if only_r:
                loop_lst = [svcomp]
            else:
                loop_lst = [svcomp, 'T']
            rayp = srad2skm(self.PArrival.ray_param)
        elif phase == 'S':
            if self.comp == 'lqt':
                loop_lst = ['L']
            else:
                loop_lst = ['Z']
            rayp = srad2skm(self.SArrival.ray_param)
        else:
            raise ValueError('Phase must be in \'P\' or \'S\'')

        if evtstr is None:
            filename = join(path, self.datestr)
        else:
            filename = join(path, evtstr)
        for comp in loop_lst:
            trrf = self.rf.select(channel='*'+comp)[0]
            header = {'evla': evla, 'evlo': evlo, 'evdp': evdp, 'mag': mag, 'baz': baz,
                      'gcarc': gcarc, 'user0': rayp, 'kuser0': 'Ray Para', 'user1': gauss, 'kuser1': 'G factor'}
            for key in kwargs:
                header[key] = kwargs[key]
            for key, value in header.items():
                trrf.stats['sac'][key] = value
            tr = SACTrace.from_obspy_trace(trrf)
            tr.b = -shift
            tr.o = 0
            tr.write(filename + '_{0}_{1}.sac'.format(phase, tr.kcmpnm[-1]))

    def judge_rf(self, shift, npts, phase='P', criterion='crust', rmsgate=None):
        if not isinstance(criterion, (str, type(None))):
            raise TypeError('criterion should be string in [\'crust\', \'mtz\']')
        elif criterion is None:
            pass
        elif criterion.lower() not in ['crust', 'mtz']:
            raise ValueError('criterion should be string in [\'crust\', \'mtz\']')
        else:
            criterion = criterion.lower()
        
        if phase == 'P' and self.comp == 'rtz':
            trrf = self.rf.select(channel='*R')[0]
        elif phase == 'P' and self.comp == 'lqt':
            trrf = self.rf.select(channel='*Q')[0]
        elif phase == 'S' and self.comp == 'lqt':
            trrf = self.rf.select(channel='*L')[0]
        elif phase == 'S' and self.comp == 'rtz':
            trrf = self.rf.select(channel='*Z')[0]        
        if trrf.stats.npts != npts:
            return False
        
        # Final RMS
        if rmsgate is not None:
            if self.method == 'iter':
                rms = self.rf[0].stats.rms[-1]
            else:
                rms = self.rf[0].stats.rms
            rmspass = rms < rmsgate
        else:
            rmspass = True
        
        # R energy
        nt1 = int(np.floor((5+shift)/trrf.stats.delta))
        reng = np.sum(np.sqrt(trrf.data[nt1:] ** 2))
        if reng < 0.1:
            rengpass = False
        else:
            rengpass = True

        # Max amplitude
        if criterion == 'crust':
            time_P1 = int(np.floor((-2+shift)/trrf.stats.delta))
            time_P2 = int(np.floor((4+shift)/trrf.stats.delta))
            max_P = np.max(trrf.data[time_P1:time_P2])
            if max_P == np.max(np.abs(trrf.data)) and max_P < 1 and rmspass and rengpass:
                return True
            else:
                return False
        elif criterion == 'mtz':
            max_deep = np.max(np.abs(trrf.data[int((30 + shift) / trrf.stats.delta):]))
            time_P1 = int(np.floor((-5 + shift) / trrf.stats.delta))
            time_P2 = int(np.floor((5 + shift) / trrf.stats.delta))
            max_P = np.max(trrf.data[time_P1:time_P2])
            if max_deep < max_P * 0.4 and rmspass and rengpass and\
                  max_P == np.max(np.abs(trrf.data)) and max_P < 1:
                return True
            else:
                return False
        elif criterion is None:
            print(rmspass, rengpass)
            return rmspass and rengpass
        else:
            pass


if __name__ == '__main__':
    pass

