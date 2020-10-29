from obspy.io.sac.sactrace import SACTrace
import numpy as np
from scipy.interpolate import interp1d, interpn
from scipy.signal import resample
from os.path import dirname, join, exists
from seispy.geo import skm2srad, sdeg2skm, rad2deg, latlon_from, \
                       asind, tand, srad2skm, km2deg
from seispy.psrayp import get_psrayp
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")


class SACStation(object):
    def __init__(self, evt_lst, only_r=False):
        """
        :param evt_lst: event list in RF data dir. Column as date string, phase, evt lat, evt lon, evt dep,
                        distance, back-azimuth, ray parameter, magnitude, gauss factor.
        """
        self.only_r = only_r
        data_path = dirname(evt_lst)
        dtype = {'names': ('evt', 'phase', 'evlat', 'evlon', 'evdep', 'dis', 'bazi', 'rayp', 'mag', 'f0'),
                 'formats': ('U20', 'U20', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
        self.event, self.phase, self.evla, self.evlo, self.evdp, self.dis, self.bazi, self.rayp, self.mag, self.f0 = \
            np.loadtxt(evt_lst, dtype=dtype, unpack=True, ndmin=1)
        # self.event = [datestr.decode() for datestr in self.event]
        # self.phase = [ph.decode() for ph in self.phase]
        self.rayp = skm2srad(self.rayp)
        self.ev_num = self.evla.shape[0]
        sample_sac = SACTrace.read(join(data_path, self.event[0] + '_' + self.phase[0] + '_R.sac'))
        self.stla = sample_sac.stla
        self.stlo = sample_sac.stlo
        self.RFlength = sample_sac.npts
        self.shift = -sample_sac.b
        self.sampling = sample_sac.delta
        self.datar = np.empty([self.ev_num, self.RFlength])
        if not only_r:
            self.datat = np.empty([self.ev_num, self.RFlength])
            for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
                sac = SACTrace.read(join(data_path, evt + '_' + ph + '_R.sac'))
                sact = SACTrace.read(join(data_path, evt + '_' + ph + '_T.sac'))
                self.datar[_i] = sac.data
                self.datat[_i] = sact.data
        else:
            for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
                sac = SACTrace.read(join(data_path, evt + '_' + ph + '_R.sac'))
                self.datar[_i] = sac.data

    def resample(self, dt):
        npts = int(self.RFlength * (self.sampling / dt)) + 1
        self.datar = resample(self.datar, npts, axis=1)
        if not self.only_r:
            self.datat = resample(self.datat, npts, axis=1)
        self.sampling = dt
        self.RFlength = npts


class DepModel(object):
    def __init__(self, YAxisRange, velmod='iasp91'):
        VelocityModel = np.loadtxt(from_file(velmod))
        self.depthsraw = VelocityModel[:, 0]
        self.vpraw = VelocityModel[:, 1]
        self.vsraw = VelocityModel[:, 2]
        self.vp = interp1d(self.depthsraw, self.vpraw)(YAxisRange)
        self.vs = interp1d(self.depthsraw, self.vsraw)(YAxisRange)
        self.depths = YAxisRange
        self.dz = np.append(0, np.diff(self.depths))
        self.R = 6371 - self.depths


def _imag2nan(arr):
    StopIndex = np.where(np.imag(arr) == 1)[0]
    if StopIndex.size != 0:
        arr[StopIndex[0]:] = np.nan
    return arr


def from_file(mode_name):
    if exists(mode_name):
        filename = mode_name
    elif exists(join(dirname(__file__), 'data', mode_name.lower()+'.vel')):
        filename = join(dirname(__file__), 'data', mode_name.lower()+'.vel')
    else:
        raise ValueError('No such velocity mode')
    return filename


def moveoutcorrect_ref(stadatar, raypref, YAxisRange, sampling, shift, velmod='iasp91'):
    """
    :param stadatar: data class of SACStation
    :param raypref: referred ray parameter in rad
    :param YAxisRange: Depth range in nd.array type
    :param sampling: dt
    :param shift: time before P
    :param velmod: Path to velocity model
    :return: Newdatar, EndIndex, x_s, x_p
    """
    dep_mod = DepModel(YAxisRange, velmod)

    x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    Tpds = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    for i in range(stadatar.ev_num):
        x_s[i] = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
        x_p[i] = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vp) ** -2)) - 1))
        Tpds[i] = np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - stadatar.rayp[i] ** 2) -
                             np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - stadatar.rayp[i] ** 2)) * (dep_mod.dz / dep_mod.R))
    Tpds_ref = np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - raypref ** 2) -
                          np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - raypref ** 2)) * (dep_mod.dz / dep_mod.R))

    Newdatar = np.zeros([stadatar.ev_num, stadatar.RFlength])
    EndIndex = np.zeros(stadatar.ev_num)
#
    for i in range(stadatar.ev_num):
        Newaxis = np.array([])
        TempTpds = Tpds[i, :]
        StopIndex = np.where(np.imag(TempTpds) == 1)[0]
        if StopIndex.size == 0:
            StopIndex = dep_mod.depths.shape[0]
        EndIndex[i] = StopIndex - 1
        Newaxis = np.append(Newaxis, np.append(np.arange(-shift, 0, sampling), 0))
        for j in np.arange(int(shift / sampling + 1), stadatar.RFlength):
            Refaxis = j * sampling - shift
            index = np.where(Refaxis <= Tpds[i, 0:StopIndex])[0]
            if index.size == 0:
                break
            Ratio = (Tpds_ref[index[0]] - Tpds_ref[index[0] - 1]) / (Tpds[i, index[0]] - Tpds[i, index[0] - 1])
            Newaxis = np.append(Newaxis, Tpds_ref[index[0] - 1] + (Refaxis - Tpds[i, index[0] - 1]) * Ratio)
        endidx = Newaxis.shape[0]
        x_new = np.arange(0, stadatar.RFlength) * sampling - shift
        Tempdata = interp1d(Newaxis, stadatar.datar[i, 0:endidx], bounds_error=False)(x_new)
        endIndice = np.where(np.isnan(Tempdata))[0]
        if endIndice.size == 0:
            New_data = Tempdata
        else:
            New_data = np.append(Tempdata[1:endIndice[0]], stadatar.datar[i, endidx+1:])
        if New_data.shape[0] < stadatar.RFlength:
            Newdatar[i] = np.append(New_data, np.zeros(stadatar.RFlength - New_data.shape[0]))
        else:
            Newdatar[i] = New_data[0: stadatar.RFlength]
    return Newdatar, EndIndex, x_s, x_p


def psrf2depth(stadatar, YAxisRange, sampling, shift, velmod='iasp91', velmod_3d=None, srayp=None):
    """
    :param stadatar:
    :param YAxisRange:
    :param sampling:
    :param shift:
    :param velmod:
    :return:
    """

    dep_mod = DepModel(YAxisRange, velmod)
    if velmod_3d is not None:
        dep_mod.vp, dep_mod.vs = interp_depth_model(velmod_3d, stadatar.stla, stadatar.stlo, YAxisRange)

    x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    Tpds = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    if srayp is None:
        for i in range(stadatar.ev_num):
            x_s[i] = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
            x_p[i] = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vp) ** -2)) - 1))
            Tpds[i] = np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - stadatar.rayp[i] ** 2) -
                                 np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - stadatar.rayp[i] ** 2)) * (dep_mod.dz / dep_mod.R))
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
            x_s[i] = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (rayp ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
            x_p[i] = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vp) ** -2)) - 1))
            Tpds[i] = np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - rayp ** 2) -
                                 np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - stadatar.rayp[i] ** 2)) * (dep_mod.dz / dep_mod.R))
    else:
        raise TypeError('srayp should be path to Ps rayp lib')

    time_axis = np.arange(0, stadatar.RFlength) * sampling - shift
    PS_RFdepth = np.zeros([stadatar.ev_num, dep_mod.depths.shape[0]])
    EndIndex = np.zeros(stadatar.ev_num)
    for i in range(stadatar.ev_num):
        TempTpds = Tpds[i, :]
        StopIndex = np.where(np.imag(TempTpds) == 1)[0]
        if StopIndex.size == 0:
            EndIndex[i] = dep_mod.depths.shape[0]
            DepthAxis = interp1d(TempTpds, dep_mod.depths, bounds_error=False)(time_axis)
        else:
            EndIndex[i] = StopIndex[0] - 1
            DepthAxis = interp1d(TempTpds[0:StopIndex], dep_mod.depths[0: StopIndex], bounds_error=False)(time_axis)

        PS_RFTempAmps = stadatar.datar[i]
        ValueIndices = np.where(np.logical_not(np.isnan(DepthAxis)))[0]

        if ValueIndices.size == 0:
            continue
        elif np.max(ValueIndices) > PS_RFTempAmps.shape[0]:
            continue
        else:
            PS_RFAmps = interp1d(DepthAxis[ValueIndices], PS_RFTempAmps[ValueIndices], bounds_error=False)(YAxisRange)
            PS_RFdepth[i] = PS_RFAmps / np.nanmax(PS_RFAmps)
    return PS_RFdepth, EndIndex, x_s, x_p


def psrf_1D_raytracing(stadatar, YAxisRange, velmod='iasp91', srayp=None):
    dep_mod = DepModel(YAxisRange, velmod)

    # x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    raylength_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplat_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    # x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    raylength_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplat_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    Tpds = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    if srayp is None:
        for i in range(stadatar.ev_num):
            x_s = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
            raylength_s[i] = (dep_mod.dz * dep_mod.R) / (np.sqrt(((dep_mod.R / dep_mod.vs) ** 2) - (stadatar.rayp[i] ** 2)) * dep_mod.vs)
            x_p = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vp) ** -2)) - 1))
            raylength_p[i] = (dep_mod.dz * dep_mod.R) / (np.sqrt(((dep_mod.R / dep_mod.vp) ** 2) - (stadatar.rayp[i] ** 2)) * dep_mod.vp)
            Tpds[i] = np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - stadatar.rayp[i] ** 2) -
                                 np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - stadatar.rayp[i] ** 2)) * (dep_mod.dz / dep_mod.R))
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
            x_s = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (rayp ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
            raylength_s[i] = (dep_mod.dz * dep_mod.R) / (np.sqrt(((dep_mod.R / dep_mod.vs) ** 2) - (rayp ** 2)) * dep_mod.vs)
            x_p = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (stadatar.rayp[i] ** 2. * (dep_mod.R / dep_mod.vp) ** -2)) - 1))
            raylength_p[i] = (dep_mod.dz * dep_mod.R) / (np.sqrt(((dep_mod.R / dep_mod.vp) ** 2) - (stadatar.rayp[i] ** 2)) * dep_mod.vp)
            Tpds[i] = np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - rayp ** 2) -
                                 np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - stadatar.rayp[i] ** 2)) * (dep_mod.dz / dep_mod.R))
            x_s = _imag2nan(x_s)
            x_p = _imag2nan(x_p)
            pplat_s[i], pplon_s[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_s))
            pplat_p[i], pplon_p[i] = latlon_from(stadatar.stla, stadatar.stlo, stadatar.bazi[i], rad2deg(x_p))
    else:
        raise TypeError('srayp should be path to Ps rayp lib')
    return pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds


def psrf_3D_raytracing(stadatar, YAxisRange, mod3d, srayp=None):
    """
Back ray trace the S wavs with a assumed ray parameter of P.

Parameters
--------------
stla: float
    The latitude of the station
stlo: float
    The longitude of the station
stadatar: object SACStation
    The data class including PRFs and more parameters
YAxisRange: array_like
    The depth array with the same intervals
mod3d: 'Mod3DPerturbation' object
    The 3D velocity model with fields of ``dep``, ``lat``,
    ``lon``, ``vp`` and ``vs``.
    """
    R = 6371 - YAxisRange
    ddepth = np.mean(np.diff(YAxisRange))
    pplat_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplat_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    pplon_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_s = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    x_p = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    Tpds = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
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
        pplat_s[i][0] = stadatar.stla
        pplon_s[i][0] = stadatar.stlo
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
        Tpds[i] = np.cumsum((np.sqrt((R / vs) ** 2 - srayps ** 2) -
                            np.sqrt((R / vp) ** 2 - stadatar.rayp[i] ** 2))
                            * (ddepth / R))
    return pplat_s, pplon_s, pplat_p, pplon_p, Tpds


def interp_depth_model(model, lat, lon, new_dep):
    #  model = np.load(modpath)
    points = [[depth, lat, lon] for depth in new_dep]
    vp = interpn((model['dep'], model['lat'], model['lon']), model['vp'], points, bounds_error=False, fill_value=None)
    vs = interpn((model['dep'], model['lat'], model['lon']), model['vs'], points, bounds_error=False, fill_value=None)
    return vp, vs


class Mod3DPerturbation:
    def __init__(self, modpath, YAxisRange, velmod='iasp91'):
        dep_mod = DepModel(YAxisRange, velmod=velmod)
        self.model = np.load(modpath)
        new1dvp = interp1d(dep_mod.depthsraw, dep_mod.vpraw)(self.model['dep'])
        new1dvs = interp1d(dep_mod.depthsraw, dep_mod.vsraw)(self.model['dep'])
        new1dvp, _, _ = np.meshgrid(new1dvp, self.model['lat'], self.model['lon'], indexing='ij')
        new1dvs, _, _ = np.meshgrid(new1dvs, self.model['lat'], self.model['lon'], indexing='ij')
        self.dvp = (self.model['vp'] - new1dvp) / new1dvp
        self.dvs = (self.model['vs'] - new1dvs) / new1dvs
        self.cvp = dep_mod.vp
        self.cvs = dep_mod.vs

    def interpdvp(self, points):
        dvp = interpn((self.model['dep'], self.model['lat'], self.model['lon']), self.dvp, points,
                      bounds_error=False, fill_value=None)
        return dvp

    def interpdvs(self, points):
        dvs = interpn((self.model['dep'], self.model['lat'], self.model['lon']), self.dvs, points,
                      bounds_error=False, fill_value=None)
        return dvs


def psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds, YAxisRange, mod3d):
    ev_num, rfdepth = raylength_p.shape
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
    time_axis = np.arange(0, stadatar.RFlength) * stadatar.sampling - stadatar.shift
    PS_RFdepth = np.zeros([stadatar.ev_num, YAxisRange.shape[0]])
    EndIndex = np.zeros(stadatar.ev_num)
    for i in range(stadatar.ev_num):
        TempTpds = Tpds[i, :]
        StopIndex = np.where(np.imag(TempTpds) == 1)[0]
        if StopIndex.size == 0:
            EndIndex[i] = YAxisRange.shape[0]
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
    stadata = SACStation('/Users/xumj/Researches/SouthYNRF/RFresult/YN001/YN001finallist.dat', only_r=True)
    vel3dmod = np.load('/Users/xumj/Researches/YN_crust/SETPvs/model_Bao_Wang.npz')
    YAxisRange = np.append(np.arange(0, 100, 1), 100)
    # PS_RFdepth, EndIndex, x_s, x_p = psrf2depth(stadata, YAxisRange, stadata.sampling, stadata.shift, velmod='iasp91', velmod_3d=vel3dmod, srayp=None)
    # plt.plot(YAxisRange, PS_RFdepth[1])
    # plt.show()
    psrf_3D_raytracing(stadata, YAxisRange, vel3dmod, srayp=None)


    '''
    # dep_mod = DepModel(YAxisRange, velmod)
    # print(dep_mod.vp.shape, dep_mod.R.shape, YAxisRange[-1])
    srayp = np.load('/Users/xumj/Researches/Tibet_MTZ/process/psrayp.npz')
    mod3d = Mod3DPerturbation('/Users/xumj/Researches/Tibet_MTZ/models/GYPSUM.npz', YAxisRange)
    pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds = psrf_1D_raytracing(stla, stlo, stadatar, YAxisRange, srayp=srayp)
    newtds = psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, Tpds, YAxisRange, mod3d)
    amp1d = time2depth(stadatar, YAxisRange, Tpds)
    amp3d = time2depth(stadatar, YAxisRange, newtds)
    plt.imshow(amp3d)
    # plt.plot(YAxisRange, amp1d[1])plt.plot(YAxisRange, amp1d[1])
    # plt.plot(YAxisRange, amp3d[1])
    # PS_RFdepth, _, x_s, x_p = psrf2depth(stadatar, YAxisRange, sampling, shift, velmod, srayp='/Users/xumj/Researches/Ps_rayp.npz')
    # mean_data = np.mean(PS_RFdepth, axis=0)

    # plt.plot(YAxisRange, mean_data)
    plt.show()
    '''

