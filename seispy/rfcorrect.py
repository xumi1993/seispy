from multiprocessing.sharedctypes import Value
from obspy.io.sac.sactrace import SACTrace
import numpy as np
from scipy.interpolate import interp1d, interpn
from scipy.signal import resample
from os.path import dirname, join, exists, isfile, abspath
from seispy.geo import skm2srad, sdeg2skm, rad2deg, latlon_from, \
                       asind, tand, srad2skm, km2deg
from seispy.psrayp import get_psrayp
from seispy.rfani import RFAni
from seispy.slantstack import SlantStack
from seispy.harmonics import Harmonics
from seispy.plotRT import plotrt as _plotrt
from seispy.plotR import plotr as _plotr
from seispy.utils import DepModel, Mod3DPerturbation
import warnings
import glob


class RFStation(object):
    def __init__(self, data_path, only_r=False, prime_comp='R'):
        """
        Class for derivative process of RFs.

        :param data_path: Path to RF data with SAC format. A finallist.dat must be in this path.
        :type data_path: str
        :param only_r: Whether reading transverse RFs, defaults to False
        :type only_r: bool, optional
        :param prime_comp: Prime component in RF filename. ``R`` or ``Q`` for PRF and ``L`` or ``Z`` for SRF
        :type prime_comp: str

        .. warning::

            This Class will be renamed to ``RFStation`` in future versions.
        """
        
        self.only_r = only_r
        self.comp = prime_comp
        self.dtype = {'names': ('event', 'phase', 'evla', 'evlo', 'evdp', 'dis', 'bazi', 'rayp', 'mag', 'f0'),
                 'formats': ('U20', 'U20', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
        if not exists(data_path):
            return
        self._chech_comp()
        if isfile(data_path):
            data_path = dirname(abspath(data_path))
        evt_lsts = glob.glob(join(data_path, '*finallist.dat'))
        if len(evt_lsts) == 0:
            raise FileNotFoundError("No such *finallist.dat in the {}".format(data_path))
        elif len(evt_lsts) > 1:
            raise ValueError("More than one finallist.dat in the {}".format(data_path))
        else:
            evt_lst = evt_lsts[0]
        self.event, self.phase, self.evla, self.evlo, self.evdp, self.dis, self.bazi, self.rayp, self.mag, self.f0 = \
            np.loadtxt(evt_lst, dtype=self.dtype, unpack=True, ndmin=1)
        self.rayp = skm2srad(self.rayp)
        self.ev_num = self.evla.shape[0]
        self.read_sample(data_path)
        self.data_prime = np.empty([self.ev_num, self.rflength])
        # self.__dict__['data{}'.format(self.comp.lower())] = np.empty([self.ev_num, self.rflength])
        # eval('self.data_prime = self.data{}'.format(self.comp.lower()))
        if not only_r:
            self.datat = np.empty([self.ev_num, self.rflength])
            for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
                sac = SACTrace.read(join(data_path, evt + '_' + ph + '_{}.sac'.format(self.comp)))
                sact = SACTrace.read(join(data_path, evt + '_' + ph + '_T.sac'))
                self.data_prime[_i] = sac.data
                self.datat[_i] = sact.data
        else:
            for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
                sac = SACTrace.read(join(data_path, evt + '_' + ph + '_{}.sac'.format(self.comp)))
                self.data_prime[_i] = sac.data
        exec('self.data{} = self.data_prime'.format(self.comp.lower()))

    def read_sample(self, data_path):
        fname = glob.glob(join(data_path, self.event[0] + '_' + self.phase[0] + '_{}.sac'.format(self.comp)))
        if len(fname) == 0:
            raise FileNotFoundError('No such files with comp of {} in {}'.format(self.comp, data_path))
        else:
            sample_sac = SACTrace.read(fname[0])
        self.staname = '{}.{}'.format(sample_sac.knetwk, sample_sac.kstnm)
        self.stla = sample_sac.stla
        self.stlo = sample_sac.stlo
        if sample_sac.stel is None:
            self.stel = 0.
        else:
            self.stel = sample_sac.stel
        self.rflength = sample_sac.npts
        self.shift = -sample_sac.b
        self.sampling = sample_sac.delta
        self.time_axis = np.arange(self.rflength) * self.sampling - self.shift

    def init_property(self, ev_num):
        self.ev_num = ev_num
        self.event = np.array(['']*ev_num)
        self.phase = np.array(['']*ev_num)
        self.evla = np.zeros(ev_num)
        self.evlo = np.zeros(ev_num)
        self.evdp = np.zeros(ev_num)
        self.mag = np.zeros(ev_num)
        self.dis = np.zeros(ev_num)
        self.bazi = np.zeros(ev_num)
        self.rayp = np.zeros(ev_num)
        self.f0 = np.zeros(ev_num)
        self.stla = 0
        self.stlo = 0
        self.stel = 0
        self.staname = ''
        self.sampling = 0
        self.rflength = 0
        # self.data_prime = np.array([])
        # exec('self.data_prime = self.data{} = np.array([])'.format(self.comp.lower()))
        # self.__dict__['data{}'.format(self.comp.lower())] = np.array([])
        self.datat = np.array([])
        # self.data_prime = eval('self.data{}'.format(self.comp.lower()))

    @classmethod
    def read_stream(cls, stream, rayp, baz, prime_comp='R'):
        """Create RFStation instance from ``obspy.Stream``

        Parameters
        ----------
        stream : ``obspy.Stream``
            _description_
        rayp : ``numpy.ndarray`` or ``float``
            Ray-parameters in s/km.
        baz : ``numpy.ndarray`` or ``float``
            Back-azimuth
        prime_comp : str, optional
             Prime component of RF, by default 'R'
        """
        if len(stream) == 0:
            raise ValueError('No such RFTrace read')
        rfsta = cls('', only_r=True, prime_comp=prime_comp)
        ev_num = len(stream)
        if isinstance(rayp, float):
            rayp = np.ones(ev_num)*rayp
        if isinstance(baz, (float, int)):
            baz = np.ones(ev_num)*baz
        if ev_num != rayp.size or ev_num != baz.size:
            raise ValueError('Array length of rayp and baz must be the same as stream')
        rfsta.init_property(ev_num)
        try:
            rfsta.staname = '{}.{}'.format(stream[0].stats.sac.knetwk, stream[0].stats.sac.kstnm)
            rfsta.stla = stream[0].stats.sac.stla
            rfsta.stlo = stream[0].stats.sac.stlo
            rfsta.stel = stream[0].stats.sac.stel
        except:
            pass
        try:
            rfsta.shift = stream[0].stats.tshift
        except:
            rfsta.shift = -stream[0].stats.sac.b
        rfsta.sampling = stream[0].stats.delta
        rfsta.rflength = stream[0].stats.npts
        rfsta.data_prime = np.zeros([rfsta.ev_num, rfsta.rflength])
        rfsta.bazi = baz
        rfsta.rayp = skm2srad(rayp)
        rfsta.time_axis = np.arange(rfsta.rflength) * rfsta.sampling - rfsta.shift
        for i, tr in enumerate(stream):
            rfsta.event[i] = '{}'.format(i)
            try:
                rfsta.f0[i] = tr.stats.f0
            except:
                rfsta.f0[i] = tr.stats.sac.user1
            rfsta.data_prime[i] = tr.data
        exec('rfsta.data{} = rfsta.data_prime'.format(rfsta.comp.lower()))
        return rfsta

    @property
    def stel(self):
        return self._stel
    
    @stel.setter
    def stel(self, value):
        if value is None:
            self._stel = 0.
        else:
            self._stel = value/1000

    def _chech_comp(self):
        if self.comp in ['R', 'Q']:
            self.prime_phase = 'P'
        elif self.comp in ['L', 'Z']:
            self.prime_phase = 'S'
        else:
            raise ValueError('prime component should be in \'R\', \'Q\', \'L\' and \'Z\'')

    def normalize(self, method='single'):
        """Normalize amplitude of each RFs.
        :param method: Method of normalization with ``single`` and ``average`` avaliable.
                     - ``single`` for normalization with max amplitude of current RF.
                     - ``average`` for normalization with average amplitude of current station.
        :type method: str, optional
        """
        if not isinstance(method, str):
            raise TypeError('\'type\' must be string, but {} type got'.format(type(method)))
        if method == 'single':
            maxamp = np.nanmax(np.abs(self.__dict__['data{}'.format(self.comp.lower())]), axis=1)
        elif method == 'average':
            amp = np.nanmax(np.abs(np.mean(self.__dict__['data{}'.format(self.comp.lower())], axis=0)))
            maxamp = np.ones(self.ev_num) * amp
        else:
            raise ValueError('\'method\' must be in \'single\' and \'average\'')
        for i in range(self.ev_num):
            self.data_prime /= maxamp[i]
            exec('self.data{} = self.data_prime'.format(self.comp.lower()))
            if not self.only_r:
                self.datat[i] /= maxamp[i]

    def resample(self, dt):
        """Resample RFs with specified dt


        Parameters
        ----------
        dt : ``float``
            Target sampling interval in sec
        """
        npts = int(self.rflength * (self.sampling / dt)) + 1
        self.data_prime = resample(self.data_prime, npts, axis=1)
        exec('self.data{} = self.data_prime'.format(self.comp.lower()))
        if not self.only_r:
            self.datat = resample(self.datat, npts, axis=1)
        self.sampling = dt
        self.rflength = npts
        self.time_axis = np.arange(npts) * dt - self.shift

    def sort(self, key='bazi'):
        """Sort RFs by keys in given ``event``, ``evla``, ``evlo``, ``evdp``,
        ``dis``, ``bazi``, ``rayp``, ``mag``, ``f0``

        :param key: key to sort, defaults to ``bazi``
        :type key: str, optional
        """
        idx = np.argsort(self.__dict__[key])
        for keyarg in self.dtype['names']:
            self.__dict__[keyarg] = self.__dict__[keyarg][idx]
        self.data_prime = self.data_prime[idx]
        exec('self.data{} = self.data_prime'.format(self.comp.lower()))
        if not self.only_r:
            self.datat = self.datat[idx]

    def moveoutcorrect(self, ref_rayp=0.06, dep_range=np.arange(0, 150), velmod='iasp91', replace=False, **kwargs):
        """Moveout correction with specified reference ray-parameter and depth

        :param ref_rayp: reference ray-parameter in s/km, defaults to 0.06
        :type ref_rayp: float, optional
        :param dep_range: Depth range used for extracting velocity in velocity model, defaults to np.arange(0, 150)
        :type dep_range: ``np.ndarray``, optional
        :param velmod: Velocity model for moveout correction. 'iasp91', 'prem' 
                      and 'ak135' is valid for internal model. Specify path to velocity model for the customized model. 
                      The format is the same as in Taup, but the depth should be monotonically increasing, defaults to 'iasp91'
        :type velmod: str, optional
        :param replace: whether replace original data, False to return new array, defaults to False
        :type replace: bool, optional

        Returns
        -------
        rf_corr: ``np.ndarray``
            Corrected RFs with component of ``RFStation.comp``

        t_corr: ``np.ndarray`` or ``None``
            Corrected RFs in transverse component. If ``only_r`` is ``True``, this variable is ``None``
        
        """
        if not self.only_r:
            t_corr, _ = moveoutcorrect_ref(self, skm2srad(ref_rayp), dep_range, chan='t', velmod=velmod, **kwargs)
        else:
            t_corr = None
        if 'datar' in self.__dict__:
            chan = 'r'
        elif 'dataz' in self.__dict__:
            chan = 'z'
        elif 'datal' in self.__dict__:
            chan = 'l'
        else:
            pass
        rf_corr, _ = moveoutcorrect_ref(self, skm2srad(ref_rayp), dep_range, chan=chan, velmod=velmod, **kwargs)
        if replace:
            self.__dict__['data{}'.format(chan)] = rf_corr
            if not self.only_r:
                self.datat = t_corr
        else:
            return rf_corr, t_corr

    def psrf2depth(self, dep_range=np.arange(0, 150), **kwargs):
        """Time-to-depth conversion with specified depth series.

        :param dep_range: Discret conversion depth, defaults to np.arange(0, 150)
        :type dep_range: :meth:`np.ndarray` , optional
        :param velmod: Velocity model for time-to-depth conversion. 'iasp91', 'prem' 
                      and 'ak135' is valid for internal model. Specify path to velocity model for the customized model. 
                      The format is the same as in Taup, but the depth should be monotonically increasing, defaults to 'iasp91'
        :type velmod: str, optional
        :param srayp: Ray-parameter lib for Ps phases, If set up to None the rayp of direct is used, defaults to None
        :type srayp: numpy.lib.npyio.NpzFile, optional
        :return: 2D array of RFs in depth
        :rtype: :meth:`np.ndarray`
        """
        self.dep_range = dep_range
        rfdepth, endindex, x_s, x_p = psrf2depth(self, dep_range, **kwargs)
        return rfdepth

    def psrf_1D_raytracing(self, dep_range=np.arange(0, 150), **kwargs):
        """1D back ray tracing to obtained Ps conversion points at discret depths

        :param dep_range: Discret conversion depth, defaults to np.arange(0, 150)
        :type dep_range: numpy.ndarray, optional
        :param velmod: Velocity model for time-to-depth conversion. 'iasp91', 'prem' 
                      and 'ak135' is valid for internal model. Specify path to velocity model for the customized model. 
                      The format is the same as in Taup, but the depth should be monotonically increasing, defaults to 'iasp91'
        :type velmod: str, optional
        :param srayp: Ray-parameter lib for Ps phases, If set up to None the rayp of direct is used, defaults to None
        :type srayp: numpy.lib.npyio.NpzFile, optional
        :return pplat_s: Latitude of conversion points
        :return pplon_s: Longitude of conversion points
        :return tps: Time difference of Ps at each depth
        :rtype: list
        """
        self.dep_range = dep_range
        pplat_s, pplon_s, _ , _, _, _, tps = psrf_1D_raytracing(self, dep_range, **kwargs)
        return pplat_s, pplon_s, tps

    def psrf_3D_raytracing(self, mod3dpath, dep_range=np.arange(0, 150), srayp=None):
        self.dep_range = dep_range
        mod3d = Mod3DPerturbation(mod3dpath, dep_range)
        pplat_s, pplon_s, _, _, tps = psrf_3D_raytracing(self, dep_range, mod3d, srayp=srayp)
        return pplat_s, pplon_s, tps

    def psrf_3D_moveoutcorrect(self, mod3dpath, **kwargs):
        warnings.warn('The fuction will be change to RFStation.psrf_3D_timecorrect in the future')
        self.psrf_3D_timecorrect(mod3dpath, **kwargs)

    def psrf_3D_timecorrect(self,  mod3dpath, dep_range=np.arange(0, 150), normalize='single', **kwargs):
        self.dep_range = dep_range
        mod3d = Mod3DPerturbation(mod3dpath, dep_range)
        pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps = psrf_1D_raytracing(self, dep_range, **kwargs)
        tps = psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps, dep_range, mod3d)
        rfdepth, _ = time2depth(self, dep_range, tps, normalize=normalize)
        return rfdepth

    def jointani(self, tb, te, tlen=3., stack_baz_val=10, rayp=0.06,
                 velmodel='iasp91', weight=[0.4, 0.4, 0.2]):
        """Eastimate crustal anisotropy with a joint method. See Liu and Niu (2012, doi: 10.1111/j.1365-246X.2011.05249.x) in detail.

        :param tb: Time before Pms for search Ps peak
        :type tb: float
        :param te: Time after Pms for search Ps peak
        :type te: float
        :param tlen: Half time length for cut out Ps phase, defaults to 3.0
        :type tlen: float, optional
        :param stack_baz_val: The interval for stacking binned by back-azimuth, defaults to 10
        :type stack_baz_val: float, optional
        :param rayp: Reference ray-parameter for moveout correction, defaults to 0.06
        :type rayp: float, optional
        :param velmodel: velocity model for moveout correction. 'iasp91', 'prem' 
                      and 'ak135' is valid for internal model. Specify path to velocity model for the customized model. 
                      The format is the same as in Taup, but the depth should be monotonically increasing, defaults to 'iasp91'
        :type velmodel: str, optional
        :param weight: Weight for three different method, defaults to [0.4, 0.4, 0.2]
        :type weight: list, optional
        :return: Dominant fast velocity direction and time delay
        :rtype: ``List``, ``List``
        """
        self.ani = RFAni(self, tb, te, tlen=tlen, rayp=rayp, model=velmodel)
        self.ani.baz_stack(val=stack_baz_val)
        best_f, best_t = self.ani.joint_ani(weight=weight)
        return best_f, best_t

    def slantstack(self, ref_dis=None, rayp_range=None, tau_range=None):
        self.slant = SlantStack(self.data_prime, self.time_axis, self.dis)
        self.slant.stack(ref_dis, rayp_range, tau_range)
        return self.slant.stack_amp

    def harmonic(self, tb=-5, te=10):
        """Harmonic decomposition for extracting anisotropic and isotropic features from the radial and transverse RFs

        :param tb: Start time relative to P, defaults to -5
        :type tb: ``float``, optional
        :param te: End time relative to P, defaults to 10
        :type te: ``float``, optional

        Returns
        -------
        harmonic_trans: numpy.ndarray, float
                Harmonic components with shape of ``(5, nsamp)``, ``nsamp = (te-tb)/RFStation.sampling``

        unmodel_trans: numpy.ndarray, float
                Unmodel components with shape same as harmonic_trans.
        """
        if self.only_r:
            raise ValueError('Transverse RFs are nessary for harmonic decomposition')
        self.harmo = Harmonics(self, tb, te)
        self.harmo.harmo_trans()
        return self.harmo.harmonic_trans, self.harmo.unmodel_trans

    def plotrt(self, outformat=None, **kwargs):
        if self.only_r:
            raise ValueError('Transverse RFs are nessary or use RFStation.plotr instead.')
        return _plotrt(self, outformat=outformat, **kwargs)

    def plotr(self, outformat=None, **kwargs):
        return _plotr(self, outformat=outformat, **kwargs)


class SACStation(RFStation):
    def __init__(self, data_path, only_r=False):
        """Class for derivative process of RFs.

        :param data_path: Path to RF data with SAC format. A finallist.dat must be in this path.
        :type data_path: str
        :param only_r: Wether only read R component, defaults to False
        :type only_r: bool, optional
        """
        super().__init__(data_path, only_r=only_r)


def _imag2nan(arr):
    StopIndex = np.where(np.imag(arr) == 1)[0]
    if StopIndex.size != 0:
        arr[StopIndex[0]:] = np.nan
    return arr


def moveoutcorrect_ref(stadatar, raypref, YAxisRange, 
                       chan='r', velmod='iasp91', sphere=True, phase=1):
    """Moveout correction refer to a specified ray-parameter
    
    :param stadatar: data class of RFStation
    :param raypref: referred ray parameter in rad
    :param YAxisRange: Depth range in nd.array type
    :param velmod: Path to velocity model
    :param chan: channel name for correction, 'r', 't'...

    :return: Newdatar, EndIndex
    """
    sampling = stadatar.sampling
    shift = stadatar.shift
    if chan == 'r':
        data = stadatar.datar
    elif chan == 't':
        data = stadatar.datat
    elif chan == 'z':
        data = stadatar.dataz
    elif chan == 'l':
        data = stadatar.datal
    else:
        raise ValueError('Field \'datar\' or \'datal\' must be in the SACStation')
    dep_mod = DepModel(YAxisRange, velmod, stadatar.stel)
    # x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    # x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    tps = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    for i in range(stadatar.ev_num):
        tps[i], _, _ = xps_tps_map(dep_mod, stadatar.rayp[i], stadatar.rayp[i], sphere=sphere, phase=phase)
    Tpds_ref, _, _ = xps_tps_map(dep_mod, raypref, raypref, sphere=sphere, phase=phase)
    Newdatar = np.zeros([stadatar.ev_num, stadatar.rflength])
    EndIndex = np.zeros(stadatar.ev_num)
#
    for i in range(stadatar.ev_num):
        Newaxis = np.array([])
        TempTpds = tps[i, :]
        StopIndex = np.where(np.imag(TempTpds) == 1)[0]
        if StopIndex.size == 0:
            StopIndex = dep_mod.depths.shape[0]
        EndIndex[i] = StopIndex - 1
        Newaxis = np.append(Newaxis, np.append(np.arange(-shift, 0, sampling), 0))
        for j in np.arange(int(shift / sampling + 1), stadatar.rflength):
            Refaxis = j * sampling - shift
            index = np.where(Refaxis <= tps[i, 0:StopIndex])[0]
            if index.size == 0:
                break
            Ratio = (Tpds_ref[index[0]] - Tpds_ref[index[0] - 1]) / (tps[i, index[0]] - tps[i, index[0] - 1])
            Newaxis = np.append(Newaxis, Tpds_ref[index[0] - 1] + (Refaxis - tps[i, index[0] - 1]) * Ratio)
        endidx = Newaxis.shape[0]
        x_new = np.arange(0, stadatar.rflength) * sampling - shift
        Tempdata = interp1d(Newaxis, data[i, 0:endidx], bounds_error=False)(x_new)
        endIndice = np.where(np.isnan(Tempdata))[0]
        if endIndice.size == 0:
            New_data = Tempdata
        else:
            New_data = np.append(Tempdata[1:endIndice[0]], data[i, endidx+1:])
        if New_data.shape[0] < stadatar.rflength:
            Newdatar[i] = np.append(New_data, np.zeros(stadatar.rflength - New_data.shape[0]))
        else:
            Newdatar[i] = New_data[0: stadatar.rflength]
    return Newdatar, EndIndex


def psrf2depth(stadatar, YAxisRange, velmod='iasp91', srayp=None, normalize='single', sphere=True, phase=1):
    """ Time-to-depth conversion with S-wave backprojection.

    :param stadatar: Data class of RFStation
    :type stadatar: :meth:`RFStation`
    :param YAxisRange: Depth range for conversion
    :type YAxisRange: ``numpy.ndarray``
    :param velmod: Velocity for conversion, whcih can be a path to velocity file, defaults to 'iasp91'
    :type velmod: str, optional
    :param srayp: ray-parameter library of conversion phases. See :meth:`seispy.psrayp` in detail, defaults to None
    :type srayp: str or :meth:`seispy.psrayp.PsRayp`, optional
    :param normalize: method of normalization, defaults to 'single'. Please refer to :meth:`RFStation.normalize`
    :type normalize: str, optional
    :param sphere: Wether do earth-flattening transformation, defaults to True
    :type sphere: bool, optional

    Returns
    ------------

    ps_rfdepth: 2-D numpy.ndarray, float
                RFs in depth with shape of ``(stadatar.ev_num, YAxisRange.size)``, ``stadatar.ev_num`` is the number of RFs in current station.
                ``YAxisRange.size`` is the size of depth axis.
    endindex: numpy.ndarray, int
                End index of each RF in depth
    x_s: 2-D numpy.ndarray, float
                Horizontal distance between station and S-wave conversion points with shape of  ``(stadatar.ev_num, YAxisRange.size)``
    x_p: 2-D numpy.ndarray, float
                Horizontal distance between station and P-wave conversion points with shape of  ``(stadatar.ev_num, YAxisRange.size)``
    """
    if exists(velmod):
        try:
            dep_mod = DepModel(YAxisRange, velmod, elevation=stadatar.stel)
        except:
            dep_mod = DepModel(YAxisRange, 'iasp91', elevation=stadatar.stel)
            try:
                velmod_3d = np.load(velmod)
                dep_mod.vp, dep_mod.vs = interp_depth_model(velmod_3d,
                                         stadatar.stla, stadatar.stlo, dep_mod.depths_elev)
            except Exception as e:
                raise FileNotFoundError('Cannot load 1D or 3D velocity model of \'{}\''.format(velmod))
    else:
        try:
            dep_mod = DepModel(YAxisRange, velmod, elevation=stadatar.stel)
        except:
            raise ValueError('Cannot recognize the velocity model of \'{}\''.format(velmod))

    x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    tps = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    if srayp is None:
        for i in range(stadatar.ev_num):
            tps[i], x_s[i], x_p[i] = xps_tps_map(dep_mod, stadatar.rayp[i], 
                stadatar.rayp[i], sphere=sphere, phase=phase)
    elif isinstance(srayp, str) or isinstance(srayp, np.lib.npyio.NpzFile):
        if isinstance(srayp, str):
            if not exists(srayp):
                raise FileNotFoundError('Ps rayp lib file was not found')
            else:
                rayp_lib = np.load(srayp)
        else:
            rayp_lib = srayp
        for i in range(stadatar.ev_num):
            rayp = get_psrayp(rayp_lib, stadatar.dis[i], stadatar.evdp[i], dep_mod.depths_elev)
            rayp = skm2srad(sdeg2skm(rayp))
            tps[i], x_s[i], x_p[i] = xps_tps_map(dep_mod,
                rayp, stadatar.rayp[i], sphere=sphere, phase=phase)
    else:
        raise TypeError('srayp should be path to Ps rayp lib')
    ps_rfdepth, endindex = time2depth(stadatar, dep_mod.depths, tps, normalize=normalize)
    return ps_rfdepth, endindex, x_s, x_p


def xps_tps_map(dep_mod, srayp, prayp, is_raylen=False, sphere=True, phase=1):
    """Calculate horizontal distance and time difference at depths

    :param dep_mod: 1D velocity model class 
    :type dep_mod: :meth:`seispy.util.DepModel`
    :param srayp: conversion phase ray-parameters
    :type srayp: float
    :param prayp: S-wave ray-parameters
    :type prayp: float
    :param is_raylen: Wether calculate ray length at depths, defaults to False
    :type is_raylen: bool, optional
    :param sphere: Wether do earth-flattening transformation, defaults to True, defaults to True
    :type sphere: bool, optional
    :param phase: Phases to calculate 1 for ``Ps``, 2 for ``PpPs``, 3 for ``PsPs+PpSs``, defaults to 1
    :type phase: int, optional

    Returns
    -----------
    If ``is_raylen = False``

    tps: 2-D numpy.ndarray, float
        RFs in depth with shape of ``(stadatar.ev_num, YAxisRange.size)``, ``stadatar.ev_num`` is the number of RFs in current station.
        ``YAxisRange.size`` is the size of depth axis.
    
    x_s: 2-D numpy.ndarray, float
        Horizontal distance between station and S-wave conversion points with shape of  ``(stadatar.ev_num, YAxisRange.size)``
    x_p: 2-D numpy.ndarray, float
        Horizontal distance between station and P-wave conversion points with shape of  ``(stadatar.ev_num, YAxisRange.size)``

    otherwise, two more variables will be returned

    raylength_s: 2-D numpy.ndarray, float
    
    raylength_p: 2-D numpy.ndarray, float   
    """
    x_s = dep_mod.radius_s(prayp, phase='S', sphere=sphere)
    x_p = dep_mod.radius_s(prayp, phase='P', sphere=sphere)
    if is_raylen:
        raylength_s = dep_mod.raylength(srayp, phase='S', sphere=sphere)
        raylength_p = dep_mod.raylength(prayp, phase='P', sphere=sphere)
    if phase == 1:
        tps = dep_mod.tpds(srayp, prayp, sphere=sphere)
    elif phase == 2:
        tps = dep_mod.tpppds(srayp, prayp, sphere=sphere)
    elif phase == 3:
        tps = dep_mod.tpspds(srayp, sphere=sphere)
    else:
        raise ValueError('Phase must be in 1 for Ps, 2 for PpPs, 3 for PsPs+PpSs')
    if dep_mod.elevation != 0:
        x_s = interp1d(dep_mod.depths_elev, x_s, bounds_error=False, fill_value=(np.nan, x_s[-1]))(dep_mod.depths)
        x_p = interp1d(dep_mod.depths_elev, x_p, bounds_error=False, fill_value=(np.nan, x_p[-1]))(dep_mod.depths)
        tps = interp1d(dep_mod.depths_elev, tps, bounds_error=False, fill_value=(np.nan, tps[-1]))(dep_mod.depths)
        if is_raylen:
            raylength_s = interp1d(dep_mod.depths_elev, raylength_s, bounds_error=False, fill_value=(np.nan, raylength_s[-1]))(dep_mod.depths)
            raylength_p = interp1d(dep_mod.depths_elev, raylength_p, bounds_error=False, fill_value=(np.nan, raylength_p[-1]))(dep_mod.depths)         
    if is_raylen:
        return tps, x_s, x_p, raylength_s, raylength_p
    else:
        return tps, x_s, x_p


def psrf_1D_raytracing(stadatar, YAxisRange, velmod='iasp91', srayp=None, sphere=True, phase=1):
    dep_mod = DepModel(YAxisRange, velmod, stadatar.stel)

    # x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    raylength_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplat_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    # x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    raylength_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplat_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    tps = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    if srayp is None:
        for i in range(stadatar.ev_num):
            tps[i], x_s, x_p, raylength_s[i], raylength_p[i] = xps_tps_map(
                dep_mod, stadatar.rayp[i], stadatar.rayp[i], is_raylen=True, sphere=sphere, phase=phase)
            pplat_s[i], pplon_s[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_s))
            pplat_p[i], pplon_p[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_p))
    elif isinstance(srayp, str) or isinstance(srayp, np.lib.npyio.NpzFile):
        if isinstance(srayp, str):
            if not exists(srayp):
                raise FileNotFoundError('Ps rayp lib file not found')
            else:
                rayp_lib = np.load(srayp)
        else:
            rayp_lib = srayp
        for i in range(stadatar.ev_num):
            rayp = get_psrayp(rayp_lib, stadatar.dis[i], stadatar.evdp[i], dep_mod.depths_elev)
            rayp = skm2srad(sdeg2skm(rayp))
            tps[i], x_s, x_p, raylength_s[i], raylength_p[i] = xps_tps_map(dep_mod, rayp, 
                                        stadatar.rayp[i], is_raylen=True, sphere=sphere, phase=phase)
            x_s = _imag2nan(x_s)
            x_p = _imag2nan(x_p)
            pplat_s[i], pplon_s[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_s))
            pplat_p[i], pplon_p[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_p))
    else:
        raise TypeError('srayp should be path to Ps rayp lib')
    return pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps


def psrf_3D_raytracing(stadatar, YAxisRange, mod3d, srayp=None, elevation=0, sphere=True):
    """
    Back ray trace the S wavs with a assumed ray parameter of P.

    :param stadatar: The data class including PRFs and more parameters
    :type stadatar: object RFStation
    :param YAxisRange: The depth array with the same intervals
    :type YAxisRange: numpy.ndarray
    :param mod3d:  The 3D velocity model with fields of ``dep``, ``lat``,
                    ``lon``, ``vp`` and ``vs``.
    :type mod3d: 'Mod3DPerturbation' object
    :param elevation: Elevation of this station relative to sea level
    :type elevation: float
    :return: pplat_s, pplon_s, pplat_p, pplon_p, tps
    :type: numpy.ndarray * 5
    """
    if sphere:
        R = 6371.0 - YAxisRange + elevation
    else:
        R = 6371.0 + elevation
    dep_range = YAxisRange.copy()
    YAxisRange -= elevation
    ddepth = np.mean(np.diff(YAxisRange))
    pplat_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplat_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    tps = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    rayps = srad2skm(stadatar.rayp)

    if isinstance(srayp, str) or isinstance(srayp, np.lib.npyio.NpzFile):
        if isinstance(srayp, str):
            if not exists(srayp):
                raise FileNotFoundError('Ps rayp lib file not found')
            else:
                rayp_lib = np.load(srayp)
        else:
            rayp_lib = srayp
    elif srayp is None:
        pass
    else:
        raise TypeError('srayp should be path to Ps rayp lib')

    for i in range(stadatar.ev_num):
        if srayp is None:
            srayps = stadatar.rayp[i]
        else:
            srayps = get_psrayp(rayp_lib, stadatar.dis[i],
                                stadatar.evdp[i], YAxisRange)
            srayps = skm2srad(sdeg2skm(srayps))
        pplat_s[i][0] = pplat_p[i][0] = stadatar.stla
        pplon_s[i][0] = pplon_p[i][0] = stadatar.stlo
        x_s[i][0] = 0
        x_p[i][0] = 0
        vs = np.zeros_like(YAxisRange)
        vp = np.zeros_like(YAxisRange)
        for j, dep in enumerate(YAxisRange[:-1]):
            vs[j] = interpn((mod3d.model['dep'], mod3d.model['lat'], mod3d.model['lon']),
                            mod3d.model['vs'], (dep, pplat_s[i, j], pplon_s[i, j]),
                            bounds_error=False, fill_value=None)
            vp[j] = interpn((mod3d.model['dep'], mod3d.model['lat'], mod3d.model['lon']),
                            mod3d.model['vp'], (dep, pplat_p[i, j], pplon_p[i, j]),
                            bounds_error=False, fill_value=None)
            x_s[i, j+1] = ddepth*tand(asind(vs[j]*rayps[i])) + x_s[i, j]
            x_p[i, j+1] = ddepth*tand(asind(vp[j]*rayps[i])) + x_p[i, j]
            pplat_s[i, j+1], pplon_s[i, j+1] = latlon_from(stadatar.stla,
                                                           stadatar.stlo,
                                                           stadatar.bazi[i],
                                                           km2deg(x_s[i, j+1]))
            pplat_p[i, j+1], pplon_p[i, j+1] = latlon_from(stadatar.stla,
                                                           stadatar.stlo,
                                                           stadatar.bazi[i],
                                                           km2deg(x_p[i, j+1]))
        tps_corr = np.cumsum((np.sqrt((R / vs) ** 2 - srayps ** 2) -
                            np.sqrt((R / vp) ** 2 - stadatar.rayp[i] ** 2))
                            * (ddepth / R))
        if elevation != 0:
            tps[i] = interp1d(YAxisRange, tps_corr)(dep_range)
    return pplat_s, pplon_s, pplat_p, pplon_p, tps


def interp_depth_model(model, lat, lon, new_dep):
    """ Interpolate Vp and Vs from 3D velocity with a specified depth range.

    Parameters
    ----------
    mod3d : :meth:`np.lib.npyio.NpzFile`
        3D velocity loaded from a ``.npz`` file
    lat : float
        Latitude of position in 3D velocity model
    lon : float
        Longitude of position in 3D velocity model
    new_dep : :meth:`np.ndarray`
        1D array of depths in km

    Returns
    -------
    Vp : :meth:`np.ndarray`
        Vp in ``new_dep``
    Vs : :meth:`np.ndarray`
        Vs in ``new_dep``
    """
    #  model = np.load(modpath)
    points = [[depth, lat, lon] for depth in new_dep]
    vp = interpn((model['dep'], model['lat'], model['lon']), model['vp'], points, bounds_error=False, fill_value=None)
    vs = interpn((model['dep'], model['lat'], model['lon']), model['vs'], points, bounds_error=False, fill_value=None)
    return vp, vs


def psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds, dep_range, mod3d):
    """ 3D time difference correction with specified ray path and 3D velocity model. 
        The input parameters can be generated with :meth:`psrf_1D_raytracing`.

    Parameters
    ----------
    pplat_s : :meth:`np.ndarray`
        2D array of latitude of S-wave in ``dep_range``, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    pplon_s : :meth:`np.ndarray`
        2D array of longitude of S-wave in ``dep_range``, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    pplat_p : :meth:`np.ndarray`
        2D array of latitude of P-wave in ``dep_range``, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    pplon_p : :meth:`np.ndarray`
        2D array of longitude of P-wave in ``dep_range``, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    raylength_s : :meth:`np.ndarray`
        2D array of ray path length of S-wave in ``dep_range``, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    raylength_p : :meth:`np.ndarray`
        2D array of ray path length of P-wave in ``dep_range``, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    Tpds : :meth:`np.ndarray`
        1D array of time difference in ``dep_range`` (``dep_range.size``)
    dep_range : :meth:`np.ndarray`
        1D array of depths in km, (``dep_range.size``)
    mod3d : :meth:`np.lib.npyio.NpzFile`
        3D velocity loaded from a ``.npz`` file

    Returns
    -------
    :meth:`np.ndarray`
        Corrected time difference in dep_range
    """
    ev_num, _ = raylength_p.shape
    timecorrections = np.zeros_like(raylength_p)
    for i in range(ev_num):
        points = np.array([dep_range, pplat_p[i], pplon_p[i]]).T
        dvp = mod3d.interpdvp(points)
        points = np.array([dep_range, pplat_s[i], pplon_s[i]]).T
        dvs = mod3d.interpdvs(points)
        dlp = raylength_p[i]
        dls = raylength_s[i]
        tmpds = (dls / (mod3d.cvs * (1 + dvs)) - dls / mod3d.cvs) - (dlp / (mod3d.cvp * (1 + dvp)) - dlp / mod3d.cvp)
        tmpds[np.isnan(tmpds)] = 0
        timecorrections[i] = np.cumsum(tmpds)
    return Tpds + timecorrections


def time2depth(stadatar, dep_range, Tpds, normalize='single'):
    """ Interpolate RF amplitude with specified time difference and depth range

    Parameters
    ----------
    stadatar : :meth:`RFStation` 
        Data class of :meth:`RFStation` 
    dep_range : :meth:`np.ndarray`
        1D array of depths in km, (``dep_range.size``)
    Tpds : :meth:`np.ndarray`
        1D array of time difference in ``dep_range`` (``dep_range.size``)
    normalize : str, optional
        Normlization option, ``'sinlge'`` and ``'average'`` are available , by default 'single'
        See :meth:`RFStation.normalize` in detail.

    Returns
    -------
    PS_RFdepth: :meth:`np.ndarray`
        2D array of RFs in depth, (:meth:`RFStation.ev_num`, ``dep_range.size``)
    end_index: :meth:`np.ndarray`
        1D array of the last digit of effective depth (:meth:`RFStation.ev_num`)
    """
    if normalize:
        stadatar.normalize(method=normalize)
    PS_RFdepth = np.zeros([stadatar.ev_num, dep_range.shape[0]])
    EndIndex = np.zeros(stadatar.ev_num).astype(int)
    for i in range(stadatar.ev_num):
        TempTpds = Tpds[i, :]
        StopIndex = np.where(np.imag(TempTpds) == 1)[0]
        if StopIndex.size == 0:
            EndIndex[i] = dep_range.size - 1
            # DepthAxis = interp1d(TempTpds, dep_range, bounds_error=False)(stadatar.time_axis)
        else:
            EndIndex[i] = StopIndex[0] - 1
            # DepthAxis = interp1d(TempTpds[0:StopIndex], dep_range[0: StopIndex], bounds_error=False)(stadatar.time_axis)

        PS_RFTempAmps = stadatar.__dict__['data{}'.format(stadatar.comp.lower())][i]
        ValueIndices = np.where(np.logical_not(np.isnan(dep_range[0:EndIndex[i]+1])))[0]
        if ValueIndices.size == 0:
            continue
        elif np.max(ValueIndices) > PS_RFTempAmps.shape[0]:
            continue
        else:
            PS_RFAmps = interp1d(stadatar.time_axis, PS_RFTempAmps, bounds_error=False)(TempTpds[0:EndIndex[i]+1])
            PS_RFdepth[i] = PS_RFAmps
    return PS_RFdepth, EndIndex


if __name__ == '__main__':
    rfsta = SACStation('/Users/xumijian/Codes/seispy-example/ex-ccp/RFresult/ZX.212')
    rfsta.jointani(2, 7, weight=[0.9, 0.1, 0.0])
    rfsta.ani.plot_polar()

