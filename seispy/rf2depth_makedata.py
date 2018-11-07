from seispy.rfcorrect import SACStation, psrf2depth, Mod3DPerturbation, psrf_1D_raytracing,\
    psrf_3D_migration, time2depth
import numpy as np
from seispy.ccppara import ccppara
from seispy.setuplog import setuplog
from seispy.geo import latlon_from, deg2km, rad2deg
from os.path import join
from scipy.io import savemat
import argparse
import sys


class Station(object):
    def __init__(self, sta_lst):
        dtype = {'names': ('station', 'evla', 'evlo'), 'formats': ('S20', 'f4', 'f4')}
        self.station, self.stla, self.stlo = np.loadtxt(sta_lst, dtype=dtype, unpack=True)
        self.station = [sta.decode() for sta in self.station]
        self.sta_num = self.stla.shape[0]


def init_mat(sta_num):
    dtype = {'names': ('Station', 'stalat', 'stalon', 'Depthrange', 'events', 'bazi', 'rayp', 'phases', 'moveout_correct',
                      'Piercelat', 'Piercelon', 'StopIndex'),
             'formats': tuple(['O']*12)}
    return np.zeros([1, sta_num], dtype=dtype)


def _convert_str_mat(instr):
    mat = np.zeros((len(instr), 1), dtype='O')
    for i in range(len(instr)):
        mat[i, 0] = np.array(np.array([instr[i]]), dtype='O')
    return mat


def makedata(cpara, use_rayp_lib=True, log=setuplog()):
    # cpara = ccppara(cfg_file)
    sta_info = Station(cpara.stalist)
    if use_rayp_lib:
        srayp = np.load(cpara.rayp_lib)
    else:
        srayp = None
    RFdepth = []
    for i in range(sta_info.stla.shape[0]):
        rfdep = {}
        evt_lst = join(cpara.rfpath, sta_info.station[i], sta_info.station[i] + 'finallist.dat')
        stadatar = SACStation(evt_lst, only_r=True)
        log.RF2depthlog.info('the {}th/{} station with {} events'.format(i + 1, sta_info.stla.shape[0], stadatar.ev_num))
        piercelat = np.zeros([stadatar.ev_num, cpara.depth_axis.shape[0]])
        piercelon = np.zeros([stadatar.ev_num, cpara.depth_axis.shape[0]])
        PS_RFdepth, end_index, x_s, x_p = psrf2depth(stadatar, cpara.depth_axis, stadatar.sampling, stadatar.shift, cpara.velmod,
                                             srayp=srayp)
        for j in range(stadatar.ev_num):
            piercelat[j], piercelon[j] = latlon_from(sta_info.stla[i], sta_info.stlo[i],
                                                     stadatar.bazi[j], rad2deg(x_s[j]))
        rfdep['Station'] = sta_info.station[i]
        rfdep['stalat'] = sta_info.stla[i]
        rfdep['stalon'] = sta_info.stlo[i]
        # rfdep['Depthrange'] = cpara.depth_axis
        # rfdep['events'] = _convert_str_mat(stadatar.event)
        rfdep['bazi'] = stadatar.bazi
        rfdep['rayp'] = stadatar.rayp
        # rfdep['phases'] = _convert_str_mat(stadatar.phase)
        rfdep['moveout_correct'] = PS_RFdepth
        rfdep['Piercelat'] = piercelat
        rfdep['Piercelon'] = piercelon
        rfdep['StopIndex'] = end_index
        RFdepth.append(rfdep)
    # savemat(cpara.depthdat, {'RFdepth': RFdepth}, oned_as='column')
    np.save(cpara.depthdat, RFdepth)


def makedata3d(cpara, velmod3d, log=setuplog()):
    mod3d = Mod3DPerturbation(velmod3d, cpara.depth_axis, velmod=cpara.velmod)
    sta_info = Station(cpara.stalist)
    srayp = np.load(cpara.rayp_lib)
    RFdepth = []
    for i in range(sta_info.stla.shape[0]):
        rfdep = {}
        evt_lst = join(cpara.rfpath, sta_info.station[i], sta_info.station[i] + 'finallist.dat')
        stadatar = SACStation(evt_lst, only_r=True)
        log.RF2depthlog.info('the {}th/{} station with {} events'.format(i + 1, sta_info.stla.shape[0], stadatar.ev_num))
        pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tpds = psrf_1D_raytracing(sta_info.stlo[i], sta_info.stla[i], stadatar,
                                                                                                cpara.depth_axis, srayp=srayp)
        newtpds = psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p,
                                    tpds, cpara.depth_axis, mod3d)
        amp3d, end_index = time2depth(stadatar, cpara.depth_axis, newtpds)
        rfdep['Station'] = sta_info.station[i]
        rfdep['stalat'] = sta_info.stla[i]
        rfdep['stalon'] = sta_info.stlo[i]
        # rfdep['Depthrange'] = cpara.depth_axis
        # rfdep['events'] = _convert_str_mat(stadatar.event)
        rfdep['bazi'] = stadatar.bazi
        rfdep['rayp'] = stadatar.rayp
        # rfdep['phases'] = _convert_str_mat(stadatar.phase)
        rfdep['moveout_correct'] = amp3d
        rfdep['Piercelat'] = pplat_s
        rfdep['Piercelon'] = pplon_s
        rfdep['StopIndex'] = end_index
        RFdepth.append(rfdep)
    savemat(cpara.depthdat, {'RFdepth': RFdepth})


def rf2depth():
    parser = argparse.ArgumentParser(description="Convert Ps RF to depth axis")
    parser.add_argument('-d', help='Path to 3d vel model in npz file', dest='vel3dpath', type=str, default='')
    parser.add_argument('cfg_file', type=str, help='Path to configure file')
    arg = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    cpara = ccppara(arg.cfg_file)
    if arg.vel3dpath == '':
        makedata(cpara)
    else:
        makedata3d(cpara, arg.vel3dpath)


if __name__ == '__main__':
    cfg_file = '/Users/xumj/Researches/Tibet_MTZ/process/paraCCP.cfg'
    cpara = ccppara(cfg_file)
    makedata(cpara)
