from seispy.rfcorrect import RFStation, psrf2depth, psrf_1D_raytracing,\
    psrf_3D_migration, time2depth, psrf_3D_raytracing
from seispy.core.pertmod import Mod3DPerturbation
import numpy as np
from seispy.ccppara import ccppara, CCPPara
from seispy.setuplog import setuplog
from seispy.geo import latlon_from, rad2deg
from os.path import join, exists
import argparse
import sys
import glob


class Station(object):
    def __init__(self, sta_lst:str):
        """
        Read station list

        :param sta_lst: Path to station list
        :type sta_lst: string
        """
        dtype = {'names': ('station', 'stla', 'stlo', 'stel'), 'formats': ('U20', 'f4', 'f4', 'f2')}
        try:
            self.station, self.stla, self.stlo, self.stel = np.loadtxt(sta_lst, dtype=dtype, unpack=True, ndmin=1)
        except:
            dtype = {'names': ('station', 'stla', 'stlo'), 'formats': ('U20', 'f4', 'f4')}
            self.station, self.stla, self.stlo = np.loadtxt(sta_lst, dtype=dtype, unpack=True, ndmin=1)
            self.stel = np.zeros(self.stla.size)
        self.sta_num = self.stla.shape[0]



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


class RFDepth():
    """Convert receiver function to depth axis
    """
    def __init__(self, cpara:CCPPara, log=setuplog(), 
                 raytracing3d=False, velmod3d=None, modfolder1d=None) -> None:
        """
        :param cpara: CCPPara object
        :type cpara: CCPPara
        :param log: Log object
        :type log: setuplog
        :param raytracing3d: If True, use 3D ray tracing to calculate the travel time
        :type raytracing3d: bool
        :param velmod3d: Path to 3D velocity model in npz file
        :type velmod3d: str
        :param modfolder1d: Folder path to 1D velocity model files with staname.vel as the file name
        :type modfolder1d: str
        """
        self.ismod1d = False
        self.cpara = cpara
        self.modfolder1d = modfolder1d
        self.log = log
        self.raytracing3d = raytracing3d
        if velmod3d is not None:
            if isinstance(velmod3d, str):
                self.mod3d = Mod3DPerturbation(velmod3d, cpara.depth_axis, velmod=cpara.velmod)
            else:
                log.RF2depthlog.error('Path to 3d velocity model should be in str')
                sys.exit(1)
        elif modfolder1d is not None:
            if isinstance(modfolder1d, str):
                if exists(modfolder1d):
                    self.ismod1d = True
                else:
                    log.RF2depthlog.error('No such folder of {}'.format(modfolder1d))
                    sys.exit(1)
            else:
                log.RF2depthlog.error('Folder to 1d velocity model files should be in str')
                sys.exit(1)
        else:
            self.ismod1d = True
        if cpara.rayp_lib is not None:
            self.srayp = np.load(cpara.rayp_lib)
        else:
            self.srayp = None
        self.sta_info = Station(cpara.stalist)
        self.rfdepth = []
        self._test_comp()

    def _test_comp(self):
        rfpath = join(self.cpara.rfpath, self.sta_info.station[0])
        self.prime_comp = ''
        for comp in ['R', 'Q', 'L', 'Z']:
            if glob.glob(join(rfpath, '*{}.sac'.format(comp))):
                self.prime_comp = comp
                break
        if not self.prime_comp:
            self.log.RF2depthlog.error('No such any RF files in \'R\',' 
                                      '\'Q\', \'L\', and \'Z\' components')
            sys.exit(1)

    def makedata(self, psphase=1):
        """Convert receiver function to depth axis

        :param psphase: 1 for Ps, 2 for PpPs, 3 for PpSs
        :type psphase: int
        """
        for i in range(self.sta_info.stla.shape[0]):
            rfpath = join(self.cpara.rfpath, self.sta_info.station[i])
            stadatar = RFStation(rfpath, only_r=True, prime_comp=self.prime_comp)
            stadatar.stel = self.sta_info.stel[i]
            stadatar.stla = self.sta_info.stla[i]
            stadatar.stlo = self.sta_info.stlo[i]
            if stadatar.prime_phase == 'P':
                sphere = True
            else:
                sphere = False
            self.log.RF2depthlog.info('the {}th/{} station with {} events'.format(i + 1, self.sta_info.stla.shape[0], stadatar.ev_num))
            if self.ismod1d:
                if self.modfolder1d is not None:
                    velmod = _load_mod(self.modfolder1d, self.sta_info.station[i])
                else:
                    velmod = self.cpara.velmod
                ps_rfdepth, end_index, x_s, _ = psrf2depth(stadatar, self.cpara.depth_axis,
                            velmod=velmod, srayp=self.cpara.rayp_lib, sphere=sphere, phase=psphase)
                piercelat, piercelon = latlon_from(self.sta_info.stla[i], self.sta_info.stlo[i],
                                                            stadatar.bazi, rad2deg(x_s))
            else:
                if self.raytracing3d:
                    pplat_s, pplon_s, pplat_p, pplon_p, newtpds = psrf_3D_raytracing(stadatar, self.cpara.depth_axis, self.mod3d, srayp=self.srayp, sphere=sphere)
                else:
                    pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps = psrf_1D_raytracing(
                        stadatar, self.cpara.depth_axis, srayp=self.srayp, sphere=sphere, phase=psphase)
                    newtpds = psrf_3D_migration(pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p,
                                                tps, self.cpara.depth_axis, self.mod3d)
                if stadatar.prime_phase == 'P':
                    piercelat, piercelon = pplat_s, pplon_s
                else:
                    piercelat, piercelon = pplat_p, pplon_p
                ps_rfdepth, end_index = time2depth(stadatar, self.cpara.depth_axis, newtpds)
            rfdep = self._write_rfdep(stadatar, ps_rfdepth, piercelat, piercelon, end_index)
            self.rfdepth.append(rfdep)
        np.save(self.cpara.depthdat, self.rfdepth)

    def _write_rfdep(self, stadata, amp, pplat, pplon, end_index):
        rfdep = {}
        rfdep['station'] = stadata.staname
        rfdep['stalat'] = stadata.stla
        rfdep['stalon'] = stadata.stlo
        rfdep['depthrange'] = self.cpara.depth_axis
        rfdep['bazi'] = stadata.bazi
        rfdep['rayp'] = stadata.rayp
        rfdep['moveout_correct'] = amp
        rfdep['piercelat'] = pplat
        rfdep['piercelon'] = pplon
        rfdep['stopindex'] = end_index
        return rfdep


def rf2depth():
    """Convert receiver function to depth axis
    """
    parser = argparse.ArgumentParser(description="Convert Ps RF to depth axis")
    parser.add_argument('-d', help='Path to 3d vel model in npz file for moveout correcting',
                        metavar='3d_velmodel_path', type=str, default='')
    parser.add_argument('-m', help='Folder path to 1d vel model files with staname.vel as the file name',
                        metavar='1d_velmodel_folder', type=str, default='')
    parser.add_argument('-r', help='Path to 3d vel model in npz file for 3D ray tracing',
                        metavar='3d_velmodel_path', type=str, default='')
    parser.add_argument('cfg_file', type=str, help='Path to configure file')
    arg = parser.parse_args()
    cpara = ccppara(arg.cfg_file)
    if arg.d != '' and arg.r != '':
        raise ValueError('Specify only 1 argument in \'-d\' and \'-r\'')
    elif arg.d != '' and arg.r == '' and arg.m == '':
        raytracing3d = False
        velmod3d = arg.d
        modfolder1d = None
    elif arg.d == '' and arg.r != '' and arg.m == '':
        raytracing3d = True
        velmod3d = arg.d
        modfolder1d = None
    elif arg.d == '' and arg.r == '' and arg.m != '':
        raytracing3d = False
        velmod3d = None
        modfolder1d = arg.m
    else:
        raytracing3d = False
        velmod3d = None
        modfolder1d = None
    rfd = RFDepth(
        cpara, raytracing3d=raytracing3d,
        velmod3d=velmod3d,
        modfolder1d=modfolder1d,
    )
    rfd.makedata()


if __name__ == '__main__':
    pass