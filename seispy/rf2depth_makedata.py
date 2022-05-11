from seispy.rfcorrect import RFStation, psrf2depth, Mod3DPerturbation, psrf_1D_raytracing,\
    psrf_3D_migration, time2depth, psrf_3D_raytracing
import numpy as np
from seispy.ccppara import ccppara
from seispy.setuplog import setuplog
from seispy.geo import latlon_from, deg2km, rad2deg
from os.path import join, dirname, exists
import argparse
import sys
import glob


class Station(object):
    def __init__(self, sta_lst):
        """
        Read station list
        """
        dtype = {'names': ('station', 'stla', 'stlo', 'stel'), 'formats': ('U20', 'f4', 'f4', 'f2')}
        try:
            self.station, self.stla, self.stlo, self.stel = np.loadtxt(sta_lst, dtype=dtype, unpack=True, ndmin=1)
        except:
            dtype = {'names': ('station', 'stla', 'stlo'), 'formats': ('U20', 'f4', 'f4')}
            self.station, self.stla, self.stlo = np.loadtxt(sta_lst, dtype=dtype, unpack=True, ndmin=1)
            self.stel = np.zeros(self.stla.size)
        self.sta_num = self.stla.shape[0]


def init_mat(sta_num):
    dtype = {'names': ('Station', 'stalat', 'stalon', 'Depthrange', 'events', 'bazi', 'rayp', 'phases', 'moveout_correct',
                      'Piercelat', 'Piercelon', 'StopIndex'),
             'formats': tuple(['O']*12)}
    return np.zeros([1, sta_num], dtype=dtype)


def _load_mod(datapath, staname):
    """Load 1D velocity model files with suffix of ".vel". The model file should be including 3 columns with depth, vp and vs.

    :param datapath: Folder name with 1D velocity model files.
    :type datapath: string
    :param staname: The station name as a part of file name of 1D velocity model files.
    :type staname: string
    """
    expresion = join(datapath, "*"+staname+"*.vel")
    modfiles = glob.glob(expresion)
    if len(modfiles) == 0:
        raise FileNotFoundError("The model file of {} were not found.".format(expresion))
    elif len(modfiles) > 1:
        raise ValueError('More then 1 file were found as the expresion: {}'.format(expresion))
    else:
        return modfiles[0]


def makedata(cpara, velmod3d=None, modfolder1d=None, log=setuplog()):
    ismod1d = False
    if velmod3d is not None:
        if isinstance(velmod3d, str):
            velmod = velmod3d
        else:
            raise ValueError('Path to 3d velocity model should be in str')
    elif modfolder1d is not None:
        if isinstance(modfolder1d, str):
            if exists(modfolder1d):
                ismod1d = True
            else:
                raise FileNotFoundError('No such folder of {}'.format(modfolder1d))
        else:
            ValueError('Path to 1d velocity model files should be in str')
    else:
        ismod1d = True

    # cpara = ccppara(cfg_file)
    sta_info = Station(cpara.stalist)
    RFdepth = []
    for i in range(sta_info.stla.shape[0]):
        rfdep = {}
        evt_lst = join(cpara.rfpath, sta_info.station[i], sta_info.station[i] + 'finallist.dat')
        stadatar = RFStation(evt_lst, only_r=True)
        stadatar.stel = sta_info.stel[i]
        stadatar.stla = sta_info.stla[i]
        stadatar.stlo = sta_info.stlo[i]
        log.RF2depthlog.info('the {}th/{} station with {} events'.format(i + 1, sta_info.stla.shape[0], stadatar.ev_num))
        piercelat = np.zeros([stadatar.ev_num, cpara.depth_axis.shape[0]])
        piercelon = np.zeros([stadatar.ev_num, cpara.depth_axis.shape[0]])
        if stadatar.prime_phase == 'P':
            sphere = True
        else:
            sphere = False
        if ismod1d:
            if modfolder1d is not None:
                velmod = _load_mod(modfolder1d, sta_info.station[i])
            else:
                velmod = cpara.velmod
        PS_RFdepth, end_index, x_s, _ = psrf2depth(stadatar, cpara.depth_axis,
                            velmod=velmod, srayp=cpara.rayp_lib, sphere=sphere, phase=cpara.phase)
        for j in range(stadatar.ev_num):
            piercelat[j], piercelon[j] = latlon_from(sta_info.stla[i], sta_info.stlo[i],
                                                     stadatar.bazi[j], rad2deg(x_s[j]))
        rfdep['station'] = sta_info.station[i]
        rfdep['stalat'] = sta_info.stla[i]
        rfdep['stalon'] = sta_info.stlo[i]
        rfdep['depthrange'] = cpara.depth_axis
        # rfdep['events'] = _convert_str_mat(stadatar.event)
        rfdep['bazi'] = stadatar.bazi
        rfdep['rayp'] = stadatar.rayp
        # rfdep['phases'] = stadatar.phase[i]
        rfdep['moveout_correct'] = PS_RFdepth
        rfdep['piercelat'] = piercelat
        rfdep['piercelon'] = piercelon
        rfdep['stopindex'] = end_index
        RFdepth.append(rfdep)
    # savemat(cpara.depthdat, {'RFdepth': RFdepth})
    np.save(cpara.depthdat, RFdepth)


def makedata3d(cpara, velmod3d, log=setuplog(), raytracing3d=True):
    mod3d = Mod3DPerturbation(velmod3d, cpara.depth_axis, velmod=cpara.velmod)
    sta_info = Station(cpara.stalist)
    if cpara.rayp_lib is not None:
        srayp = np.load(cpara.rayp_lib)
    else:
        srayp = None
    RFdepth = []
    for i in range(sta_info.stla.shape[0]):
        rfdep = {}
        evt_lst = join(cpara.rfpath, sta_info.station[i], sta_info.station[i] + 'finallist.dat')
        stadatar = RFStation(evt_lst, only_r=True)
        stadatar.stel = sta_info.stel[i]
        stadatar.stla = sta_info.stla[i]
        stadatar.stlo = sta_info.stlo[i]
        if stadatar.prime_phase == 'P':
            sphere = True
        else:
            sphere = False
        log.RF2depthlog.info('the {}th/{} station with {} events'.format(i + 1, sta_info.stla.shape[0], stadatar.ev_num))
        if raytracing3d:
            pplat_s, pplon_s, pplat_p, pplon_p, newtpds = psrf_3D_raytracing(stadatar, cpara.depth_axis, mod3d, srayp=srayp, sphere=sphere)
        else:
            pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps = psrf_1D_raytracing(
                stadatar, cpara.depth_axis, srayp=srayp, sphere=sphere, phase=cpara.phase)
            newtpds = psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p,
                                        tps, cpara.depth_axis, mod3d)
        amp3d, end_index = time2depth(stadatar, cpara.depth_axis, newtpds)
        rfdep['station'] = sta_info.station[i]
        rfdep['stalat'] = sta_info.stla[i]
        rfdep['stalon'] = sta_info.stlo[i]
        rfdep['depthrange'] = cpara.depth_axis
        # rfdep['events'] = _convert_str_mat(stadatar.event)
        rfdep['bazi'] = stadatar.bazi
        rfdep['rayp'] = stadatar.rayp
        # rfdep['phases'] = _convert_str_mat(stadatar.phase)
        rfdep['moveout_correct'] = amp3d
        rfdep['piercelat'] = pplat_s
        rfdep['piercelon'] = pplon_s
        rfdep['stopindex'] = end_index
        RFdepth.append(rfdep)
    np.save(cpara.depthdat, RFdepth)


def rf2depth():
    parser = argparse.ArgumentParser(description="Convert Ps RF to depth axis")
    parser.add_argument('-d', help='Path to 3d vel model in npz file for moveout correcting',
                        metavar='3d_velmodel_path', type=str, default='')
    parser.add_argument('-m', help='Folder path to 1d vel model files with staname.vel as the file name',
                        metavar='1d_velmodel_folder', type=str, default='')
    parser.add_argument('-r', help='Path to 3d vel model in npz file for 3D ray tracing',
                        metavar='3d_velmodel_path', type=str, default='')
    parser.add_argument('cfg_file', type=str, help='Path to configure file')
    arg = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    cpara = ccppara(arg.cfg_file)
    if arg.d != '' and arg.r != '':
        raise ValueError('Specify only 1 argument in \'-d\' and \'-r\'')
    elif arg.d != '' and arg.r == '' and arg.m == '':
        makedata3d(cpara, arg.d, raytracing3d=False)
    elif arg.d == '' and arg.r != '' and arg.m == '':
        makedata3d(cpara, arg.r, raytracing3d=True)
    elif arg.d == '' and arg.r == '' and arg.m != '':
        makedata(cpara, modfolder1d=arg.m)
    else:
        makedata(cpara)


if __name__ == '__main__':
    cfg_file = '/Users/xumj/Researches/Tibet_MTZ/process/paraCCP.cfg'
    vel3d_file = '/Users/xumj/Researches/Tibet_MTZ/models/GYPSUM.npz'
    cpara = ccppara(cfg_file)
    makedata3d(cpara, vel3d_file)
