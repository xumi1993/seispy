from obspy.io.sac.sactrace import SACTrace
import numpy as np
from scipy.interpolate import interp1d, interpn
from scipy.signal import resample
from os.path import dirname, join, exists, basename, isfile, abspath
from seispy.geo import skm2srad, sdeg2skm, rad2deg, latlon_from, \
                       asind, tand, srad2skm, km2deg
from seispy.psrayp import get_psrayp
from seispy.rfani import RFAni
from seispy.slantstack import SlantStack
# import matplotlib.pyplot as plt
from seispy.utils import DepModel, tpds, Mod3DPerturbation
import warnings
import glob


class SACStation(object):
    def __init__(self, data_path, only_r=False):
        """
        Class for derivative process of RFs.

        :param data_path: Path to RF data with SAC format. A finallist.dat must be in this path.
        :type data_path: str
        :param only_r: [description], defaults to False
        :type only_r: bool, optional

        .. warning::

            This Class will be renamed to ``RFStation`` in future versions.
        """
        
        self.only_r = only_r
        if isfile(data_path):
            data_path = dirname(abspath(data_path))
        self.staname = basename(abspath(data_path))
        evt_lsts = glob.glob(join(data_path, '*finallist.dat'))
        if len(evt_lsts) == 0:
            raise FileNotFoundError("No such *finallist.dat in the {}".format(data_path))
        elif len(evt_lsts) > 1:
            raise ValueError("More than one finallist.dat in the {}".format(data_path))
        else:
            evt_lst = evt_lsts[0]
        dtype = {'names': ('evt', 'phase', 'evlat', 'evlon', 'evdep', 'dis', 'bazi', 'rayp', 'mag', 'f0'),
                 'formats': ('U20', 'U20', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
        self.event, self.phase, self.evla, self.evlo, self.evdp, self.dis, self.bazi, self.rayp, self.mag, self.f0 = \
            np.loadtxt(evt_lst, dtype=dtype, unpack=True, ndmin=1)
        # self.event = [datestr.decode() for datestr in self.event]
        # self.phase = [ph.decode() for ph in self.phase]
        self.rayp = skm2srad(self.rayp)
        self.ev_num = self.evla.shape[0]
        try:
            sample_sac = SACTrace.read(join(data_path, self.event[0] + '_' + self.phase[0] + '_R.sac'))
            self.comp = 'R'
        except:
            sample_sac = SACTrace.read(join(data_path, self.event[0] + '_' + self.phase[0] + '_Q.sac'))
            self.comp = 'Q'
        self.stla = sample_sac.stla
        self.stlo = sample_sac.stlo
        self.stel = sample_sac.stel
        self.rflength = sample_sac.npts
        self.shift = -sample_sac.b
        self.sampling = sample_sac.delta
        self.time_axis = np.arange(self.rflength) * self.sampling - self.shift
        self.datar = np.empty([self.ev_num, self.rflength])
        if not only_r:
            self.datat = np.empty([self.ev_num, self.rflength])
            for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
                sac = SACTrace.read(join(data_path, evt + '_' + ph + '_{}.sac'.format(self.comp)))
                sact = SACTrace.read(join(data_path, evt + '_' + ph + '_T.sac'))
                self.datar[_i] = sac.data
                self.datat[_i] = sact.data
        else:
            for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
                sac = SACTrace.read(join(data_path, evt + '_' + ph + '_{}.sac'.format(self.comp)))
                self.datar[_i] = sac.data

    @property
    def stel(self):
        return self._stel
    
    @stel.setter
    def stel(self, value):
        if value is None:
            self._stel = 0.
        else:
            self._stel = value/1000

    def normalize(self):
        maxamp = np.nanmax(np.abs(self.datar), axis=1)
        for i in range(self.ev_num):
            self.datar[i] /= maxamp[i]
            if not self.only_r:
                self.datat[i] /= maxamp[i]

    def resample(self, dt):
        npts = int(self.rflength * (self.sampling / dt)) + 1
        self.datar = resample(self.datar, npts, axis=1)
        if not self.only_r:
            self.datat = resample(self.datat, npts, axis=1)
        self.sampling = dt
        self.rflength = npts
        self.time_axis = np.arange(npts) * dt - self.shift

    def moveoutcorrect(self, ref_rayp=0.06, dep_range=np.arange(0, 150), velmod='iasp91', replace=False):
        if not self.only_r:
            t_corr, _ = moveoutcorrect_ref(self, skm2srad(ref_rayp), dep_range, chan='t', velmod=velmod)
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
        rf_corr, _ = moveoutcorrect_ref(self, skm2srad(ref_rayp), dep_range, chan=chan, velmod=velmod)
        if replace:
            self.__dict__['data{}'.format(chan)] = rf_corr
            if not self.only_r:
                self.datat = t_corr
        else:
            return rf_corr, t_corr

    def psrf2depth(self, dep_range=np.arange(0, 150), velmod='iasp91', srayp=None):
        self.dep_range = dep_range
        rfdepth, _, _, _ = psrf2depth(self, dep_range, velmod=velmod, srayp=srayp)
        return rfdepth

    def psrf_1D_raytracing(self, dep_range=np.arange(0, 150), velmod='iasp91', srayp=None):
        self.dep_range = dep_range
        pplat_s, pplon_s, _ , _, _, _, tpds = psrf_1D_raytracing(self, dep_range, velmod=velmod, srayp=srayp)
        return pplat_s, pplon_s, tpds

    def psrf_3D_raytracing(self, mod3dpath, dep_range=np.arange(0, 150), srayp=None):
        self.dep_range = dep_range
        mod3d = Mod3DPerturbation(mod3dpath, dep_range)
        pplat_s, pplon_s, _, _, tpds = psrf_3D_raytracing(self, dep_range, mod3d, srayp=srayp)
        return pplat_s, pplon_s, tpds

    def psrf_3D_moveoutcorrect(self, mod3dpath, **kwargs):
        warnings.warn('The fuction will be change to RFStation.psrf_3D_timecorrect in the future')
        self.psrf_3D_timecorrect(mod3dpath, **kwargs)

    def psrf_3D_timecorrect(self,  mod3dpath, dep_range=np.arange(0, 150), velmod='iasp91', srayp=None):
        self.dep_range = dep_range
        mod3d = Mod3DPerturbation(mod3dpath, dep_range)
        pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds = psrf_1D_raytracing(self, dep_range, velmod=velmod, srayp=srayp)
        tps = psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds, dep_range, mod3d)
        rfdepth, _ = time2depth(self, dep_range, tps)
        return rfdepth

    def jointani(self, tb, te, tlen=3., stack_baz_val=10, rayp=0.06,
                 model='iasp91', weight=[0.4, 0.4, 0.2]):
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
        :param model: velocity model for moveout correction. 'iasp91', 'prem' 
                      and 'ak135' is valid for internal model. Specify path to velocity model for the customized model. 
                      The format is the same as in Taup, but the depth should be monotonically increasing, defaults to 'iasp91'
        :type model: str, optional
        :param weight: Weight for three different method, defaults to [0.4, 0.4, 0.2]
        :type weight: list, optional
        :return: Dominant fast velocity direction and time delay
        :rtype: list, list
        """
        self.ani = RFAni(self, tb, te, tlen=tlen, rayp=rayp, model=model)
        self.ani.baz_stack(val=stack_baz_val)
        best_f, best_t = self.ani.joint_ani(weight=weight)
        return best_f, best_t

    def slantstack(self, ref_dis=None, rayp_range=None, tau_range=None):
        self.slant = SlantStack(self.datar, self.time_axis, self.dis)
        self.slant.stack(ref_dis, rayp_range, tau_range)
        return self.slant.stack_amp


class SRFStation(SACStation):
    def __init__(self, path):
        fnames = sorted(glob.glob(join(path, '*S_*.sac')))
        self.ev_num = len(fnames)
        self.rayp = np.zeros(self.ev_num)
        self.bazi = np.zeros(self.ev_num)
        self.dis = np.zeros(self.ev_num)
        self.evla = np.zeros(self.ev_num)
        self.evlo = np.zeros(self.ev_num)
        self.evdp = np.zeros(self.ev_num)
        sample_sac = SACTrace.read(fnames[0])
        self.sampling = sample_sac.delta
        self.stla = sample_sac.stla
        self.stlo = sample_sac.stlo
        self.shift = -sample_sac.b
        self.rflength = sample_sac.npts
        self.time_axis = np.arange(self.rflength) * self.sampling - self.shift
        self.datal = np.zeros([self.ev_num, self.rflength])
        for i, sacfile in enumerate(fnames):
            sac = SACTrace.read(sacfile)
            self.rayp[i] = sac.user0
            self.bazi[i] = sac.baz
            self.dis[i] = sac.gcarc
            self.evla[i] = sac.evla
            self.evlo[i] = sac.evlo
            self.evdp[i] = sac.evdp
            self.datal[i] = sac.data
        self.rayp = skm2srad(self.rayp)

class RFStation(SACStation):
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
                       chan='r', velmod='iasp91'):
    """Moveout correction refer to a specified ray-parameter
    
    :param stadatar: data class of SACStation
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
        tps[i], _, _ = xps_tps_map(dep_mod, stadatar.rayp[i], stadatar.rayp[i])
    Tpds_ref, _, _ = xps_tps_map(dep_mod, raypref, raypref)
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


def psrf2depth(stadatar, YAxisRange, velmod='iasp91', srayp=None):
    """
    :param stadatar:
    :param YAxisRange:
    :param sampling:
    :param shift:
    :param velmod:
    :return:
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
            tps[i], x_s[i], x_p[i] = xps_tps_map(dep_mod, stadatar.rayp[i], stadatar.rayp[i])
    elif isinstance(srayp, str) or isinstance(srayp, np.lib.npyio.NpzFile):
        if isinstance(srayp, str):
            if not exists(srayp):
                raise FileNotFoundError('Ps rayp lib file was not found')
            else:
                rayp_lib = np.load(srayp)
        else:
            rayp_lib = srayp
        for i in range(stadatar.ev_num):
            rayp = get_psrayp(rayp_lib, stadatar.dis[i], stadatar.evdp[i], dep_mod.depths)
            rayp = skm2srad(sdeg2skm(rayp))
            tps[i], x_s[i], x_p[i] = xps_tps_map(dep_mod, rayp, stadatar.rayp[i])
    else:
        raise TypeError('srayp should be path to Ps rayp lib')
    ps_rfdepth, endindex = time2depth(stadatar, dep_mod.depths, tps)
    return ps_rfdepth, endindex, x_s, x_p

    # time_axis = np.arange(0, stadatar.rflength) * sampling - shift
    # PS_RFdepth = np.zeros([stadatar.ev_num, dep_mod.depths.shape[0]])
    # EndIndex = np.zeros(stadatar.ev_num)
    # for i in range(stadatar.ev_num):
    #     TempTpds = Tpds[i, :]
    #     StopIndex = np.where(np.imag(TempTpds) == 1)[0]
    #     if StopIndex.size == 0:
    #         EndIndex[i] = dep_mod.depths.shape[0]
    #         DepthAxis = interp1d(TempTpds, dep_mod.depths, bounds_error=False)(time_axis)
    #     else:
    #         EndIndex[i] = StopIndex[0] - 1
    #         DepthAxis = interp1d(TempTpds[0:StopIndex], dep_mod.depths[0: StopIndex], bounds_error=False)(time_axis)

    #     PS_RFTempAmps = stadatar.datar[i]
    #     ValueIndices = np.where(np.logical_not(np.isnan(DepthAxis)))[0]

    #     if ValueIndices.size == 0:
    #         continue
    #     elif np.max(ValueIndices) > PS_RFTempAmps.shape[0]:
    #         continue
    #     else:
    #         PS_RFAmps = interp1d(DepthAxis[ValueIndices], PS_RFTempAmps[ValueIndices], bounds_error=False)(YAxisRange)
    #         PS_RFdepth[i] = PS_RFAmps / np.nanmax(PS_RFAmps)
    # return PS_RFdepth, EndIndex, x_s, x_p


def xps_tps_map(dep_mod, srayp, prayp, is_raylen=False):
    x_s = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (srayp ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
    x_p = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (prayp ** 2. * (dep_mod.R / dep_mod.vp) ** -2)) - 1))
    if is_raylen:
        raylength_s = (dep_mod.dz * dep_mod.R) / (np.sqrt(((dep_mod.R / dep_mod.vs) ** 2) - (srayp ** 2)) * dep_mod.vs)
        raylength_p = (dep_mod.dz * dep_mod.R) / (np.sqrt(((dep_mod.R / dep_mod.vs) ** 2) - (prayp ** 2)) * dep_mod.vs)
    tps = tpds(dep_mod, srayp, prayp)
    if dep_mod.elevation != 0:
        x_s = interp1d(dep_mod.depths_elev, x_s)(dep_mod.depths)
        x_p = interp1d(dep_mod.depths_elev, x_p)(dep_mod.depths)
        tps = interp1d(dep_mod.depths_elev, tps)(dep_mod.depths)
        if is_raylen:
            raylength_s = interp1d(dep_mod.depths_elev, raylength_s)(dep_mod.depths)
            raylength_p = interp1d(dep_mod.depths_elev, raylength_p)(dep_mod.depths)         
    if is_raylen:
        return tps, x_s, x_p, raylength_s, raylength_p
    else:
        return tps, x_s, x_p


def psrf_1D_raytracing(stadatar, YAxisRange, velmod='iasp91', srayp=None):
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
            tps[i], x_s, x_p, raylength_s[i], raylength_p[i] = xps_tps_map(dep_mod, stadatar.rayp[i], stadatar.rayp[i], is_raylen=True)
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
            rayp = get_psrayp(rayp_lib, stadatar.dis[i], stadatar.evdp[i], dep_mod.depths)
            rayp = skm2srad(sdeg2skm(rayp))
            tps[i], x_s, x_p, raylength_s[i], raylength_p[i] = xps_tps_map(dep_mod, rayp, stadatar.rayp[i], is_raylen=True)
            x_s = _imag2nan(x_s)
            x_p = _imag2nan(x_p)
            pplat_s[i], pplon_s[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_s))
            pplat_p[i], pplon_p[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_p))
    else:
        raise TypeError('srayp should be path to Ps rayp lib')
    return pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps


def psrf_3D_raytracing(stadatar, YAxisRange, mod3d, srayp=None, elevation=0):
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
    R = 6371.0 - YAxisRange + elevation
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
    #  model = np.load(modpath)
    points = [[depth, lat, lon] for depth in new_dep]
    vp = interpn((model['dep'], model['lat'], model['lon']), model['vp'], points, bounds_error=False, fill_value=None)
    vs = interpn((model['dep'], model['lat'], model['lon']), model['vs'], points, bounds_error=False, fill_value=None)
    return vp, vs


def psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds, YAxisRange, mod3d):
    ev_num, _ = raylength_p.shape
    timecorrections = np.zeros_like(raylength_p)
    for i in range(ev_num):
        points = np.array([YAxisRange, pplat_p[i], pplon_p[i]]).T
        dvp = mod3d.interpdvp(points)
        points = np.array([YAxisRange, pplat_s[i], pplon_s[i]]).T
        dvs = mod3d.interpdvs(points)
        dlp = raylength_p[i]
        dls = raylength_s[i]
        tmpds = (dls / (mod3d.cvs * (1 + dvs)) - dls / mod3d.cvs) - (dlp / (mod3d.cvp * (1 + dvp)) - dlp / mod3d.cvp)
        tmpds[np.isnan(tmpds)] = 0
        timecorrections[i] = np.cumsum(tmpds)
    return Tpds + timecorrections


def time2depth(stadatar, YAxisRange, Tpds):
    time_axis = np.arange(0, stadatar.rflength) * stadatar.sampling - stadatar.shift
    PS_RFdepth = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    EndIndex = np.zeros(stadatar.ev_num)
    for i in range(stadatar.ev_num):
        TempTpds = Tpds[i, :]
        StopIndex = np.where(np.imag(TempTpds) == 1)[0]
        if StopIndex.size == 0:
            EndIndex[i] = YAxisRange.shape[0] - 1
            DepthAxis = interp1d(TempTpds, YAxisRange, bounds_error=False)(time_axis)
        else:
            EndIndex[i] = StopIndex[0] - 1
            DepthAxis = interp1d(TempTpds[0:StopIndex], YAxisRange[0: StopIndex], bounds_error=False)(time_axis)

        PS_RFTempAmps = stadatar.datar[i]
        ValueIndices = np.where(np.logical_not(np.isnan(DepthAxis)))[0]

        if ValueIndices.size == 0:
            continue
        elif np.max(ValueIndices) > PS_RFTempAmps.shape[0]:
            continue
        else:
            PS_RFAmps = interp1d(DepthAxis[ValueIndices], PS_RFTempAmps[ValueIndices], bounds_error=False)(YAxisRange)
            PS_RFdepth[i] = PS_RFAmps / np.nanmax(PS_RFAmps)
    return PS_RFdepth, EndIndex


if __name__ == '__main__':
    rfsta = SACStation('/Users/xumijian/Codes/seispy-example/ex-ccp/RFresult/ZX.212/ZX.212finallist.dat')
    rfsta.jointani(2, 7, weight=[0.9, 0.1, 0.0])
    rfsta.ani.plot_polar()

