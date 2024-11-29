"""
class RF2depth : process cal RF2depth
class sta_part, sta_full, _RFInd

"""
from os.path import join, exists
import argparse
import sys
import glob
from collections import namedtuple
from logging import Logger

import numpy as np

from seispy.core.depmodel import _load_mod, DepModel
from seispy.rfcorrect import RFStation, psrf2depth, psrf_1D_raytracing,\
    psrf_3D_migration, time2depth, psrf_3D_raytracing
from seispy.core.pertmod import Mod3DPerturbation
import numpy as np
from seispy.ccppara import ccppara, CCPPara
from seispy.setuplog import SetupLog
from seispy.geo import latlon_from, rad2deg

# sta_part is not used
sta_part = namedtuple('sta_part', ['station', 'stla', 'stlo'])
sta_full = namedtuple('sta_full', ['station', 'stla', 'stlo', 'stel'])

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

    def __getitem__(self,index):
        """
        allow for sta[index]
        """
        if hasattr(self,'stel'):
            return sta_full(self.station[index],
                            self.stla[index], self.stlo[index], self.stel[index])
        else:
            return sta_full(self.station[index],
                            self.stla[index], self.stlo[index], 0)

    def __iter__(self):
        """
        allow for sta in stalist
        """
        for index in range(self.stla.shape[0]):
            if hasattr(self, 'stel'):
                yield sta_full(self.station[index],
                               self.stla[index], self.stlo[index], self.stel[index])
            else:
                yield sta_full(self.station[index], self.stla[index], self.stlo[index], 0)
    def __len__(self):
        return self.station.shape[0]



class RFDepth():
    """Convert receiver function to depth axis
    """
    def __init__(self, cpara:CCPPara, log:Logger=SetupLog().RF2depthlog,
                 raytracing3d=False, velmod3d=None, modfolder1d=None) -> None:
        """
        Convert receiver function to depth axis
        
        :param cpara: CCPPara object
        :type cpara: CCPPara
        :param log: Log object
        :type log: logging.Logger
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
                log.error('Path to 3d velocity model should be in str')
                sys.exit(1)
        elif modfolder1d is not None:
            if isinstance(modfolder1d, str):
                if exists(modfolder1d):
                    self.ismod1d = True
                else:
                    log.error('No such folder of {}'.format(modfolder1d))
                    sys.exit(1)
            else:
                log.error('Folder to 1d velocity model files should be in str')
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
            self.log.error('No such any RF files in \'R\',' 
                                      '\'Q\', \'L\', and \'Z\' components')
            sys.exit(1)

    def makedata(self, psphase=1):
        """Convert receiver function to depth axis

        :param psphase: 1 for Ps, 2 for PpPs, 3 for PpSs
        :type psphase: int
        """
        for _i, _sta in enumerate(self.sta_info):
            rfpath = join(self.cpara.rfpath, _sta.station)
            stadatar = RFStation(rfpath, only_r=True, prime_comp=self.prime_comp)
            stadatar.stel = _sta.stel
            stadatar.stla = _sta.stla
            stadatar.stlo = _sta.stlo
            if stadatar.prime_phase == 'P':
                sphere = True
            else:
                sphere = False
            self.log.info('the {}th/{} station with {} events'.format(_i + 1, len(self.sta_info), stadatar.ev_num))

            #### 1d model for each station
            if self.ismod1d:
                if self.modfolder1d is not None:
                    velmod = _load_mod(self.modfolder1d, _sta.station)
                else:
                    velmod = self.cpara.velmod

                ps_rfdepth, end_index, x_s, _ = psrf2depth(stadatar, self.cpara.depth_axis,
                            velmod=velmod, srayp=self.cpara.rayp_lib,
                           sphere=sphere, phase=psphase)

                piercelat, piercelon = np.zeros_like(x_s, dtype=np.float64), np.zeros_like(x_s, dtype=np.float64)

                for j in range(stadatar.ev_num):
                    piercelat[j], piercelon[j] = latlon_from(_sta.stla, _sta.stlo,
                                                            stadatar.bazi[j], rad2deg(x_s[j]))
            else:
                ### 3d model interp
                if self.raytracing3d:
                    pplat_s, pplon_s, pplat_p, pplon_p, newtpds = psrf_3D_raytracing(
                        stadatar, self.cpara.depth_axis, self.mod3d, srayp=self.srayp, sphere=sphere
                    )
                else:
                    pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p, tps = psrf_1D_raytracing(
                        stadatar, self.cpara.depth_axis, velmod=self.cpara.velmod, srayp=self.srayp, sphere=sphere, phase=psphase
                    )
                    newtpds = psrf_3D_migration(
                        pplat_s, pplon_s, pplat_p, pplon_p, raylength_s, raylength_p,
                        tps, self.cpara.depth_axis, self.mod3d
                    )
                if stadatar.prime_phase == 'P':
                    piercelat, piercelon = pplat_s, pplon_s
                else:
                    piercelat, piercelon = pplat_p, pplon_p

                ps_rfdepth, end_index = time2depth(stadatar, self.cpara.depth_axis, newtpds)

            self._write_rfdep(stadatar, ps_rfdepth, piercelat, piercelon, end_index)
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
        self.rfdepth.append(rfdep)


def rf2depth():
    """
    CLI for Convert receiver function to depth axis
    There's  4 branch provided to do RF 2 depth conversion

    1. only -d :do moveout correction
    2. only -r : do raytracing but no moveout correction
    3. -d and -r : do moveout correction and raytracing
    4. -m : use {staname}.vel file for RF2depth conversion

    """
    parser = argparse.ArgumentParser(description="Convert Ps RF to depth axis")
    parser.add_argument('-d',
                        help='Path to 3d vel model in npz file for moveout correcting',
                        metavar='3d_velmodel_path', type=str, default='')
    parser.add_argument('-m',
                        help='Folder path to 1d vel model files with staname.vel as the file name',
                        metavar='1d_velmodel_folder', type=str, default='')
    parser.add_argument('-r',
                        help='Path to 3d vel model in npz file for 3D ray tracing',
                        metavar='3d_velmodel_path', type=str, default='')
    parser.add_argument('cfg_file',
                        help='Path to configure file',
                        metavar='ccp.cfg', type=str)
    arg = parser.parse_args()
    cpara = ccppara(arg.cfg_file)
    if arg.d != '' and arg.r != '':
        #### print help issue
        raise ValueError('Specify only 1 argument in \'-d\' and \'-r\'')
    elif arg.d != '' and arg.r == '' and arg.m == '':
        #### do 3d moveout correction but 1d rf2depth
        raytracing3d = False
        velmod3d = arg.d
        modfolder1d = None
    elif arg.d == '' and arg.r != '' and arg.m == '':
        #### do 3d raytraying
        raytracing3d = True
        velmod3d = arg.d
        modfolder1d = None
    elif arg.d == '' and arg.r == '' and arg.m != '':
        #### use multiple 1d velmod for time2depth convertion
        raytracing3d = False
        velmod3d = None
        modfolder1d = arg.m
    else:
        ### all last, use default 1D vel model do time2depth conversion
        raytracing3d = False
        velmod3d = None
        modfolder1d = None
    rfd = RFDepth(
        cpara,
        raytracing3d=raytracing3d,
        velmod3d=velmod3d,
        modfolder1d=modfolder1d,
    )
    rfd.makedata()


if __name__ == '__main__':
    pass