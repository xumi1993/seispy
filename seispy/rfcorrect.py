from obspy.io.sac.sactrace import SACTrace
import numpy as np
from scipy.interpolate import interp1d
from os.path import dirname, join, exists
from seispy.geo import skm2srad, sdeg2skm
from seispy.psrayp import get_psrayp
import matplotlib.pyplot as plt


class SACStation(object):
    def __init__(self, evt_lst):
        """
        :param evt_lst: event list in RF data dir. Column as date string, phase, evt lat, evt lon, evt dep,
                        distance, back-azimuth, ray parameter, magnitude, gauss factor.
        """
        data_path = dirname(evt_lst)
        dtype = {'names': ('evt', 'phase', 'evlat', 'evlon', 'evdep', 'dis', 'bazi', 'rayp', 'mag', 'f0'),
                 'formats': ('U20', 'U20', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
        self.event, self.phase, self.evla, self.evlo, self.evdp, self.dis, self.bazi, self.rayp, self.mag, self.f0 = \
            np.loadtxt(evt_lst, dtype=dtype, unpack=True)
        # self.event = [datestr.decode() for datestr in self.event]
        # self.phase = [ph.decode() for ph in self.phase]
        self.rayp = skm2srad(self.rayp)
        self.ev_num = len(self.event)
        sample_sac = SACTrace.read(join(data_path, self.event[0] + '_' + self.phase[0] + '_R.sac'))
        self.RFlength = sample_sac.npts
        self.shift = -sample_sac.b
        self.sampling = sample_sac.delta
        self.datar = np.empty([self.ev_num, self.RFlength])
        self.datat = np.empty([self.ev_num, self.RFlength])
        for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
            sac = SACTrace.read(join(data_path, evt + '_' + ph + '_R.sac'))
            sact = SACTrace.read(join(data_path, evt + '_' + ph + '_T.sac'))
            self.datar[_i] = sac.data
            self.datat[_i] = sact.data


class DepModel(object):
    def __init__(self, YAxisRange, velmod):
        VelocityModel = np.loadtxt(velmod)
        Depths = VelocityModel[:, 0]
        Vp = VelocityModel[:, 1]
        Vs = VelocityModel[:, 2]
        self.vp = interp1d(Depths, Vp)(YAxisRange)
        self.vs = interp1d(Depths, Vs)(YAxisRange)
        self.depths = YAxisRange
        self.dz = np.append(0, np.diff(self.depths))
        self.R = 6371 - self.depths


def moveoutcorrect_ref(stadatar, raypref, YAxisRange, sampling, shift, velmod):
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


def psrf2depth(stadatar, YAxisRange, sampling, shift, velmod, srayp=None):
    """
    :param stadatar:
    :param YAxisRange:
    :param sampling:
    :param shift:
    :param velmod:
    :return:
    """
    dep_mod = DepModel(YAxisRange, velmod)

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


if __name__ == '__main__':
    velmod = '/Users/xumj/Researches/YN_crust/gradient_study/iasp91.vel'
    RFlength = 13001
    sampling = 0.01
    shift = 10
    raypref = skm2srad(0.065)
    lst = '/Volumes/xumj3/YNRF/RFresult/MC17/MC17finallist.dat'
    YAxisRange = np.append(np.arange(0, 200, 0.5), 200)
    # dep_mod = DepModel(YAxisRange, velmod)
    # print(dep_mod.vp.shape, dep_mod.R.shape, YAxisRange[-1])

    stadatar = SACStation(lst)
    PS_RFdepth, _, x_s, x_p = psrf2depth(stadatar, YAxisRange, sampling, shift, velmod, srayp='/Users/xumj/Researches/Ps_rayp.npz')
    print(PS_RFdepth.shape)
    mean_data = np.mean(PS_RFdepth, axis=0)

    plt.plot(YAxisRange, mean_data)
    plt.show()

