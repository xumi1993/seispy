import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.rotate import rotate2zne, rotate_zne_lqt
from scipy.signal import resample
from os.path import join
from seispy.decon import RFTrace
from seispy.geo import snr, srad2skm, rotateSeisENtoTR, \
                       rssq, extrema
from obspy.signal.trigger import recursive_sta_lta
import glob


class NotEnoughComponent(Exception):
    def __init__(self, matchkey):
        self.matchkey = matchkey
    def __str__(self):
        return '{}'.format(self.matchkey)

class TooMoreComponents(Exception):
    def __init__(self, matchkey):
        self.matchkey = matchkey
    def __str__(self):
        return '{}'.format(self.matchkey)

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


class EQ(object):
    def __init__(self, pathname, datestr, suffix='SAC'):
        """Class for processing event data with 3 components, which read SAC files of ``pathname*datastr*suffix`` 

        :param pathname: Directory to SAC files
        :type pathname: string
        :param datestr: date part in filename, e.g., ``2021.122.12.23.40``
        :type datestr: string
        :param suffix: suffix for SAC files, defaults to 'SAC'
        :type suffix: str, optional
        """
        self.datestr = datestr
        self.filestr = join(pathname, '*' + datestr + '*' + suffix)
        if glob.glob(self.filestr):  
            self.st = obspy.read(self.filestr)
            self._check_comp()
            self.st.sort()
            self.set_comp()
        else:
            self.st = obspy.Stream()
        self.rf = obspy.Stream()
        self.timeoffset = 0
        self.rms = np.array([0])
        self.it = 0
        self.trigger_shift = 0
        self.inc_correction = 0
    
    def _check_comp(self):
        if len(self.st) < 3:
            channel = ' '.join([tr.stats.channel for tr in self.st])
            raise NotEnoughComponent('Sismogram must be in 3 components, but there are only channel {} of {}'.format(channel, self.datestr))
        elif len(self.st) > 3:
            raise TooMoreComponents('{} has more than 3 components, please select to delete redundant seismic components'.format(self.datestr))
        else:
            pass

    def readstream(self):
        self.rf = obspy.Stream()
        self.st = obspy.read(self.filestr)
        self._check_comp()
        self.st.sort()
        self.set_comp()

    def cleanstream(self):
        self.st = None
        self.rf = None

    @classmethod
    def from_stream(cls, stream):
        eq = cls('', '')
        eq.st = stream
        eq.datestr = eq.st[0].stats.starttime.strftime('%Y.%j.%H.%M.%S')
        eq._check_comp()
        eq.st.sort()
        eq.set_comp()
        return eq

    def set_comp(self):
        if self.st.select(channel='*[E2]'):
            self.comp = 'enz'
        elif self.st.select(channel='*R'):
            self.comp = 'rtz'
        elif self.st.select(channel='*Q'):
            self.comp = 'lqt'
        else:
            raise ValueError('No such component in R, E or Q')

    def channel_correct(self, switchEN=False, reverseE=False, reverseN=False):
        """_summary_

        :param switchEN: _description_, defaults to False
        :type switchEN: bool, optional
        :param reverseE: _description_, defaults to False
        :type reverseE: bool, optional
        :param reverseN: _description_, defaults to False
        :type reverseN: bool, optional
        """
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

    def write(self, path, evt_datetime):
        for tr in self.st:
            sac = SACTrace.from_obspy_trace(tr)
            sac.b = 0
            sac.o = evt_datetime - tr.stats.starttime
            fname = join(path, '{}.{}.{}.{}.SAC'.format(
                tr.stats.network, tr.stats.station,
                tr.stats.starttime.strftime('%Y.%j.%H%M%S'),
                tr.stats.channel
            ))
            sac.write(fname)

    def detrend(self):
        self.st.detrend(type='linear')
        self.st.detrend(type='constant')
        self.fix_channel_name()

    def filter(self, freqmin=0.05, freqmax=1, order=4):
        self.st.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=order, zerophase=True)

    def get_arrival(self, model, evdp, dis, phase='P'):
        arrivals = model.get_travel_times(evdp, dis, phase_list=[phase])
        if not arrivals:
            raise ValueError('The phase of {} is not exists'.format(phase))
        if len(arrivals) > 1:
            raise ValueError('More than one phase were calculated with distance of {} and focal depth of {}'.format(dis, evdp))
        else:
            # self.arrival = arrivals[0]
            self.arr_time = arrivals[0].time
            self.rayp = arrivals[0].ray_param
            self.inc = arrivals[0].incident_angle
            self.phase = phase

    def search_inc(self, bazi):
        inc_range = np.arange(0.1, 90, 0.1)
        s_range = self.trim(20, 20, isreturn=True)
        power = np.zeros(inc_range.shape[0])
        for i in range(len(inc_range)):
            l_comp, _, _ = rotate_zne_lqt(s_range[2].data, s_range[1].data, s_range[0].data, bazi, inc_range[i])
            power[i] = np.mean(l_comp ** 2)

        real_inc_idx = np.argmin(power)
        real_inc = inc_range[real_inc_idx]
        self.inc_correction = real_inc - self.inc
        self.inc = real_inc

    def search_baz(self, bazi, time_b=10, time_e=20, offset=90):
        p_arr = self.arr_correct(write_to_sac=False)
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

    def rotate(self, baz, inc=None, method='NE->RT', search_inc=False, baz_shift=0):
        bazi = np.mod(baz + baz_shift, 360)
        if inc is None:
            if self.phase[-1] == 'S' and search_inc:
                inc = self.search_inc(bazi)
        else:
            self.inc = inc

        if method == 'NE->RT':
            self.comp = 'rtz'
            self.st.rotate('NE->RT', back_azimuth=bazi)
        elif method == 'RT->NE':
            self.st.rotate('RT->NE', back_azimuth=bazi)
        elif method == 'ZNE->LQT':
            self.comp = 'lqt'
            self.st.rotate('ZNE->LQT', back_azimuth=bazi, inclination=self.inc)
        elif method == 'LQT->ZNE':
            self.st.rotate('LQT->ZNE', back_azimuth=bazi, inclination=self.inc)
        else:
            pass

    def snr(self, length=50):
        st_noise = self.trim(length, 0, isreturn=True)
        st_signal = self.trim(0, length, isreturn=True)
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
    
    def get_time_offset(self, event_time=None):
        if event_time is not None and not isinstance(event_time, obspy.core.utcdatetime.UTCDateTime):
            raise TypeError('Event time should be UTCDateTime type in obspy')
        elif event_time is None:
            self.timeoffset = self.st[2].stats.sac.b - self.st[2].stats.sac.o
        else:
            self.timeoffset = self.st[2].stats.starttime - event_time

    def arr_correct(self, write_to_sac=True):
        """
        offset = sac.b - real o
        """
        correct_time = self.arr_time - self.timeoffset
        if write_to_sac:
            for tr in self.st:
                tr.stats.sac.t0 = correct_time
                tr.stats.sac.kt0 = self.phase

        return correct_time

    def _get_time(self, time_before, time_after):
        arr = self.arr_correct(write_to_sac=False)
        t1 = self.st[2].stats.starttime + (arr + self.trigger_shift - time_before)
        t2 = self.st[2].stats.starttime + (arr + self.trigger_shift + time_after)
        return t1, t2

    def phase_trigger(self, time_before, time_after, stl=5, ltl=10):
        t1, t2 = self._get_time(time_before, time_after)
        self.st_pick = self.st.copy().trim(t1, t2)
        if len(self.st_pick) == 0:
            return
        if self.phase[-1] == 'P':
            tr = self.st_pick.select(channel='*Z')[0]
        else:
            tr = self.st_pick.select(channel='*T')[0]
        df = tr.stats.sampling_rate
        cft = recursive_sta_lta(tr.data, int(stl*df), int(ltl*df))
        n_trigger = np.argmax(np.diff(cft)[int(ltl*df):])+int(ltl*df)
        self.t_trigger = t1 + n_trigger/df
        self.trigger_shift = n_trigger/df - time_before

    def trim(self, time_before, time_after, isreturn=False):
        """
        offset = sac.b - real o
        """
        t1, t2 = self._get_time(time_before, time_after)
        if isreturn:
            return self.st.copy().trim(t1, t2)
        else:
            self.st.trim(t1, t2)

    def deconvolute(self, shift, time_after, f0=2.0, method='iter', only_r=False,
                    itmax=400, minderr=0.001, wlevel=0.05, target_dt=None):
        """Deconvolution

        Parameters
        ----------
        shift : float
            Time shift before P arrival
        time_after : float
            Time length after P arrival
        f0 : float or list, optional
            Gaussian factors, by default 2.0
        method : str, optional
            method for deconvolution in ``iter`` or ``water``, by default ``iter``
        only_r : bool, optional
            Whether only calculate RF in prime component, by default False
        itmax : int, optional
            Maximum iterative number, valid for method of ``iter``, by default 400
        minderr : float, optional
            Minium residual error, valid for method of ``iter``, by default 0.001
        wlevel : float, optional
            Water level, valid for method of ``water``, by default 0.05
        target_dt : None or float, optional
            Time delta for resampling, by default None
        """
        self.method = method
        if isinstance(f0, (int, float)):
            f0 = [f0]
        for ff in f0:
            if method == 'iter':
                kwargs = {'method': method,
                        'f0': ff,
                        'tshift': shift,
                        'itmax': itmax,
                        'minderr': minderr}
            elif method == 'water':
                kwargs = {'method': method,
                        'f0': ff,
                        'tshift': shift,
                        'wlevel': wlevel}
            else:
                raise ValueError('method must be in \'iter\' or \'water\'')

            if self.phase[-1] == 'P':
                self.decon_p(**kwargs)
                if not only_r:
                    self.decon_p(tcomp=True, **kwargs)
            else:
                # TODO: if 'Q' not in self.rf[1].stats.channel or 'L' not in self.rf[2].stats.channel:
                #     raise ValueError('Please rotate component to \'LQT\'')
                self.decon_s(**kwargs)
            if target_dt is not None:
                for tr in self.rf:
                    if tr.stats.delta != target_dt:
                        tr.data = resample(tr.data, int((shift + time_after)/target_dt+1))
                        tr.stats.delta = target_dt
    
    def decon_p(self, tshift, tcomp=False, **kwargs):
        if self.comp == 'lqt':
            win = self.st.select(channel='*L')[0]
            if tcomp:
                uin = self.st.select(channel='*T')[0]
            else:
                uin = self.st.select(channel='*Q')[0]
                uin.data *= -1
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
        # win.data[0:int((tshift-4)/win.stats.delta)] = 0
        uout = RFTrace.deconvolute(uin, win, phase='S', tshift=tshift, **kwargs)
        uout.data = np.flip(uout.data)
        self.rf.append(uout)

    def saverf(self, path, evtstr=None, shift=0, evla=-12345., evlo=-12345., evdp=-12345., mag=-12345.,
               gauss=0, baz=-12345., gcarc=-12345., only_r=False, **kwargs):
        if self.phase[-1] == 'P':
            if self.comp == 'lqt':
                svcomp = 'Q'
            else:
                svcomp = 'R'
            if only_r:
                loop_lst = [svcomp]
            else:
                loop_lst = [svcomp, 'T']
            rayp = srad2skm(self.rayp)
        elif self.phase[-1] == 'S':
            if self.comp == 'lqt':
                loop_lst = ['L']
            else:
                loop_lst = ['Z']
            rayp = srad2skm(self.rayp)
        else:
            pass
        if evtstr is None:
            filename = join(path, self.datestr)
        else:
            filename = join(path, evtstr)
        for comp in loop_lst:
            trrfs = self.rf.select(channel='*'+comp)
            try:
                trrf = [tr for tr in trrfs if tr.stats.f0 == gauss][0]
            except:
                ValueError('No such gauss factor of {} in calculated RFs'.format(gauss))
            header = {'evla': evla, 'evlo': evlo, 'evdp': evdp, 'mag': mag, 'baz': baz,
                      'gcarc': gcarc, 'user0': rayp, 'kuser0': 'Ray Para', 'user1': gauss, 'kuser1': 'G factor'}
            for key in kwargs:
                header[key] = kwargs[key]
            for key, value in header.items():
                trrf.stats['sac'][key] = value
            tr = SACTrace.from_obspy_trace(trrf)
            tr.b = -shift
            tr.a = 0
            tr.ka = self.phase
            tr.write(filename + '_{0}_{1}.sac'.format(self.phase, tr.kcmpnm[-1]))

    def s_condition(self, trrf, shift):
        nt0 = int(np.floor((shift)/trrf.stats.delta))
        nt25 = int(np.floor((shift+25)/trrf.stats.delta))
        if rssq(trrf.data[nt0:nt25]) > rssq(trrf.data[nt25:]):
            return True
        else:
            return False

    def judge_rf(self, gauss, shift, npts, criterion='crust', rmsgate=None):
        if self.phase[-1] == 'P' and self.comp == 'rtz':
            trrfs = self.rf.select(channel='*R')
        elif self.phase[-1] == 'P' and self.comp == 'lqt':
            trrfs = self.rf.select(channel='*Q')
        elif self.phase[-1] == 'S' and self.comp == 'lqt':
            trrfs = self.rf.select(channel='*L')
        elif self.phase[-1] == 'S' and self.comp == 'rtz':
            trrfs = self.rf.select(channel='*Z')
        for tr in trrfs:      
            if tr.stats.npts != npts:
                return False
        try:
            trrf = [tr for tr in trrfs if tr.stats.f0 == gauss][0]
        except:
            ValueError('No such gauss factor of {} in calculated RFs'.format(gauss))
        
        # All points are NaN
        if np.isnan(trrf.data).all():
            return False
        
        if np.isinf(trrf.data).any():
            return False
        
        # Final RMS
        if rmsgate is not None:
            if self.method == 'iter':
                rms = trrf.stats.rms[-1]
            else:
                rms = trrf.stats.rms
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
                return True and rengpass
            else:
                return False
        elif criterion == 'mtz':
            max_deep = np.max(np.abs(trrf.data[int((30 + shift) / trrf.stats.delta):]))
            time_P1 = int(np.floor((-5 + shift) / trrf.stats.delta))
            time_P2 = int(np.floor((5 + shift) / trrf.stats.delta))
            max_P = np.max(trrf.data[time_P1:time_P2])
            if max_deep < max_P * 0.4 and rmspass and rengpass and\
                  max_P == np.max(np.abs(trrf.data)) and max_P < 1:
                return True and rengpass
            else:
                return False
        elif criterion == "lab":
            return self.s_condition(trrf, shift) and rengpass
        elif criterion is None:
            return rmspass and rengpass
        else:
            pass


if __name__ == '__main__':
    pass

