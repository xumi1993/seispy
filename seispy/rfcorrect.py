import obspy
from obspy.io.sac.sactrace import SACTrace
import numpy as np
from scipy.interpolate import interp1d
from os.path import dirname, join
from seispy.geo import skm2srad


class SACStation(object):
    def __init__(self, evt_lst):
        """
        :param evt_lst: event list in RF data dir. Column as date string, phase, evt lat, evt lon, evt dep,
                        distance, back-azimuth, ray parameter, magnitude, gauss factor.
        """
        data_path = dirname(evt_lst)
        dtype = {'names': ('evt', 'phase', 'evlat', 'evlon', 'evdep', 'dis', 'bazi', 'rayp', 'mag', 'f0'),
                 'formats': ('S20', 'S10', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')}
        self.event, self.phase, self.evla, self.evlo, self.evdp, self.dis, self.bazi, self.rayp, self.mag, self.f0 = \
            np.loadtxt(evt_lst, dtype=dtype, unpack=True)
        self.event = [datestr.decode() for datestr in self.event]
        self.phase = [ph.decode() for ph in self.phase]
        self.rayp = skm2srad(self.rayp)
        self.ev_num = len(self.event)
        self.RFlength = SACTrace.read(join(data_path, self.event[0] + '_' + self.phase[0] + '_R.sac')).npts
        self.data = np.empty([self.ev_num, self.RFlength])
        for _i, evt, ph in zip(range(self.ev_num), self.event, self.phase):
            sac = SACTrace.read(join(data_path, evt + '_' + ph + '_R.sac'))
            self.data[_i] = sac.data


class DepModel(object):
    def __init__(self, YAxisRange, velmod):
        VelocityModel = np.loadtxt(velmod)
        Depths = VelocityModel[:, 0]
        Vp = VelocityModel[:, 2]
        Vs = VelocityModel[:, 1]
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
        New_data = np.array([])
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
        Tempdata = interp1d(Newaxis, stadatar.data[i, 0:endidx], bounds_error=False)(x_new)
        endIndice = np.where(np.isnan(Tempdata))[0]
        if endIndice.size == 0:
            New_data = Tempdata
        else:
            New_data = np.append(Tempdata[1:endIndice[0]], stadatar.data[i, endidx+1:])
        if New_data.shape[0] < stadatar.RFlength:
            Newdatar[i] = np.append(New_data, np.zeros(stadatar.RFlength - New_data.shape[0]))
        else:
            Newdatar[i] = New_data[0: stadatar.RFlength]
    return Newdatar, EndIndex, x_s, x_p


if __name__ == '__main__':
    velmod = '/Users/xumj/Researches/YN_crust/gradient_study/IASP91.vel'
    RFlength = 13001
    sampling = 0.01
    shift = 10
    raypref = skm2srad(0.065)
    lst = '/Volumes/xumj2/CXRF/RFresult/PZH01/T1.PZH01finallist.dat'
    YAxisRange = np.append(np.arange(0, 300, 0.5), 300)
    # dep_mod = DepModel(YAxisRange, velmod)
    # print(dep_mod.vp.shape, dep_mod.R.shape, YAxisRange[-1])

    stadatar = SACStation(lst)
    Newdatar, _, _, _ = moveoutcorrect_ref(stadatar, raypref, YAxisRange, sampling, shift, velmod)
    mean_data = np.mean(Newdatar, axis=0)

