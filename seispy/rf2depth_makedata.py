from seispy.rfcorrect import SACStation, psrf2depth
import numpy as np
from seispy.ccppara import ccppara
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


def makedata(cfg_file):
    cpara = ccppara(cfg_file)
    sta_info = Station(cpara.stalist)
    RFdepth = init_mat(sta_info.stla.shape[0])
    for i in range(sta_info.stla.shape[0]):
        print('the {}th station----------------'.format(i+1))
        evt_lst = join(cpara.rfpath, sta_info.station[i], sta_info.station[i] + 'finallist.dat')
        stadatar = SACStation(evt_lst)
        piercelat = np.zeros([stadatar.ev_num, cpara.depth_axis.shape[0]])
        piercelon = np.zeros([stadatar.ev_num, cpara.depth_axis.shape[0]])
        PS_RFdepth, end_index, x_s, x_p = psrf2depth(stadatar, cpara.depth_axis, stadatar.sampling, stadatar.shift, cpara.velmod,
                                             srayp=cpara.rayp_lib)
        for j in range(stadatar.ev_num):
            piercelat[j], piercelon[j] = latlon_from(sta_info.stla[i], sta_info.stlo[i],
                                                     stadatar.bazi[j], deg2km(rad2deg(x_s[j])))
        RFdepth[0, i]['Station'] = sta_info.station[i]
        RFdepth[0, i]['stalat'] = sta_info.stla[i]
        RFdepth[0, i]['stalon'] = sta_info.stlo[i]
        RFdepth[0, i]['Depthrange'] = cpara.depth_axis
        RFdepth[0, i]['events'] = _convert_str_mat(stadatar.event)
        RFdepth[0, i]['bazi'] = stadatar.bazi
        RFdepth[0, i]['rayp'] = stadatar.rayp
        RFdepth[0, i]['phases'] = _convert_str_mat(stadatar.phase)
        RFdepth[0, i]['moveout_correct'] = PS_RFdepth.T
        RFdepth[0, i]['Piercelat'] = piercelat.T
        RFdepth[0, i]['Piercelon'] = piercelon.T
        RFdepth[0, i]['StopIndex'] = end_index
    savemat(cpara.depthdat, {'RFdepth': RFdepth}, oned_as='column')


def rf2depth():
    parser = argparse.ArgumentParser(description="Convert Ps RF to depth axis")
    parser.add_argument('cfg_file', type=str, help='Path to configure file')
    arg = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    makedata(arg.cfg_file)


if __name__ == '__main__':
    cfg_file = '/Users/xumj/Codes/seispy/Scripts/paraCCP.cfg'
    makedata(cfg_file)