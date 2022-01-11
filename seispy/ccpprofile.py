import numpy as np
import seispy
from seispy.geo import km2deg, deg2km, latlon_from, geoproject, sind, rad2deg, skm2srad
from seispy.setuplog import setuplog
from seispy.distaz import distaz
from seispy.rfcorrect import DepModel
from seispy.rf2depth_makedata import Station
from seispy.ccppara import ccppara, CCPPara
from scikits.bootstrap import ci
from seispy.ccp3d import boot_bin_stack
from seispy.utils import check_stack_val, read_rfdep, create_center_bin_profile
from os.path import exists, dirname, basename, join


def line_proj(lat1, lon1, lat2, lon2):
    daz = distaz(lat1, lon1, lat2, lon2)
    az1_begin = (daz.baz - 90) % 360
    az2_begin = (daz.baz + 90) % 360
    az1_end = (daz.az - 90) % 360
    az2_end = (daz.az + 90) % 360
    return az1_begin, az2_begin, az1_end, az2_end, daz


def fix_filename(filename, typ='dat'):
    dname = dirname(filename)
    if not exists(dname) and dname != '':
        raise FileExistsError('internal error')
    bname = basename(filename)
    sp_name = bname.split('.')
    if len(sp_name) == 1:
        return filename + '.' + typ
    if sp_name[-1] != typ:
        nname = '.'.join(sp_name[0:-1]) + '.' + typ
        return join(dirname(filename), nname)
    else:
        return filename


def bin_shape(cpara):
    if cpara.bin_radius is None:
        depmod = DepModel(cpara.stack_range)
        fzone = km2deg(np.sqrt(0.5*cpara.domperiod*depmod.vs*cpara.stack_range))
    else:
        fzone = np.ones_like(cpara.stack_range) * km2deg(cpara.bin_radius)
    return fzone


def init_profile(lat1, lon1, lat2, lon2, val):
    """ Initial bins along a profile with given position of two points.

    :param lat1: The latitude of the start point
    :type lat1: float
    :param lon1: The lontitude of the start point
    :type lon1: float
    :param lat2: The latitude of the end point
    :type lat2: float
    :param lon2: The lontitude of the end point
    :type lon2: float
    :param val: The interval between two points in km
    :type val: float
    :return: The location of bins (bin_loca), and length between each bin and the start point (profile_range)

    The bin_loca is positions of bins with a numpy.array with two column. The profile_range is distance between bin center and the start point with an 1D numpy.array.
    :rtype: (numpy.array, numpy.array)
    """
    azi = distaz(lat1, lon1, lat2, lon2).baz
    dis = distaz(lat1, lon1, lat2, lon2).delta
    profile_range = np.arange(0, deg2km(dis), val)
    lat_loca, lon_loca = latlon_from(lat1, lon1, azi, km2deg(profile_range))
    bin_loca = np.zeros([lat_loca.shape[0], 2])
    bin_loca = np.vstack((lat_loca, lon_loca)).T
    return bin_loca, profile_range


class CCPProfile():
    def __init__(self, cfg_file=None, log=None):
        if log is None:
            self.logger = setuplog()
        else:
            self.logger = log
        if cfg_file is None:
            self.cpara = CCPPara()
        elif isinstance(cfg_file, str):
            self.load_para(cfg_file)
        else:
            raise ValueError('cfg_file must be str format.')
        self.stack_data = []

    def load_para(self, cfg_file):
        try:
            self.cpara = ccppara(cfg_file)
        except Exception as e:
            self.logger.CCPlog.error('Cannot open configure file {}'.format(cfg_file))
            raise FileNotFoundError('{}'.format(e))
        try:
            self.stack_mul = check_stack_val(self.cpara.stack_val, self.cpara.dep_val)
        except Exception as e:
            self.logger.CCPlog.error('{}'.format(e))
            raise ValueError('{}'.format(e))

    def read_rfdep(self):
        """Read RFdepth file

        :raises FileNotFoundError: Not Found RFdepth file
        """
        self.logger.CCPlog.info('Loading RFdepth data from {}'.format(self.cpara.depthdat))
        try:
            self.rfdep = read_rfdep(self.cpara.depthdat)
        except FileNotFoundError as e:
            self.logger.CCPlog.error('{}'.format(e))
            raise FileNotFoundError('Cannot open file of {}'.format(self.cpara.depthdat))
    
    def initial_profile(self):
        """Initialize bins of profile
        """
        self.read_rfdep()
        if exists(self.cpara.stack_sta_list):
            self.stations = Station(self.cpara.stack_sta_list)
        if self.cpara.adaptive:
            if not exists(self.cpara.stack_sta_list):
                raise ValueError('Adaptive binning depends on existence of {} order by stations'.format(
                    self.cpara.stack_sta_list))
            lat, lon, self.profile_range = create_center_bin_profile(self.stations, self.cpara.slide_val)
            self.bin_loca = np.vstack((lat, lon)).T
        else:
            self.bin_loca, self.profile_range = init_profile(*self.cpara.line, self.cpara.slide_val)
        self.fzone = bin_shape(self.cpara)
        self._get_sta()
        self._select_sta()
    
    def _get_sta(self):
        """Read station info from rfdep
        """
        self.staname = [sta['station'] for sta in self.rfdep]
        self.stalst = np.array([[sta['stalat'], sta['stalon']] for sta in self.rfdep])

    def _select_sta(self):
        if exists(self.cpara.stack_sta_list):
            self.logger.CCPlog.info('Use stacking stations in {}'.format(self.cpara.stack_sta_list))
            self.idxs = []
            for sta in self.stations.station:
                try:
                    idx = self.staname.index(sta)
                    self.idxs.append(idx)
                except ValueError:
                    self.logger.CCPlog.warning('{} does not in RFdepth structure'.format(sta))
                if self.cpara.shape == 'rect' and self.cpara.adaptive == False:
                    self._pierce_project(self.rfdep[idx])
        elif self.cpara.width is None and self.cpara.shape == 'circle':
            dep_mod = DepModel(self.cpara.stack_range)
            x_s = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (skm2srad(0.085) ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
            dis = self.fzone[-1] + rad2deg(x_s[-1]) + 0.3
            # self.idxs = self._proj_sta(dis)
            self.idxs = []
            for i, bin_info in enumerate(self.bin_loca):
                self.idxs.append(np.where(distaz(bin_info[0], bin_info[1], self.stalst[:, 0], self.stalst[:, 1]).delta <= dis)[0])
        elif self.cpara.width is not None and self.cpara.shape == 'rect':
            self.logger.CCPlog.info('Select stations within {} km perpendicular to the profile'.format(self.cpara.width))
            self.idxs = self._proj_sta(self.cpara.width)
            self._write_sta()  
        else:
            if self.cpara.width is not None:
                self.logger.CCPlog.error('Width of profile was set, the bin shape of {} is invalid'.format(self.cpara.shape))
                raise ValueError('Width of profile was set, the bin shape of {} is invalid'.format(self.cpara.shape))
            else:
                self.logger.CCPlog.error('Width of profile was not set, the bin shape of {} is invalid'.format(self.cpara.shape))
                raise ValueError('Width of profile was not set, the bin shape of {} is invalid'.format(self.cpara.shape))

    def _write_sta(self):
        with open(self.cpara.stack_sta_list, 'w') as f:
            for idx in self.idxs:
                f.write('{}\t{:.3f}\t{:.3f}\n'.format(self.staname[idx],
                self.stalst[idx,0], self.stalst[idx, 1]))
                self._pierce_project(self.rfdep[idx])      

    def _proj_sta(self, width):
        az1_begin, az2_begin, az1_end, az2_end, daz = line_proj(*self.cpara.line)
        az_sta_begin = distaz(self.stalst[:, 0], self.stalst[:, 1], self.cpara.line[0], self.cpara.line[1]).az
        az_sta_end = distaz(self.stalst[:, 0], self.stalst[:, 1], self.cpara.line[2], self.cpara.line[3]).az
        if 0 <= daz.baz < 90 or 270 <= daz.baz < 360:
            tmp_idx_begin = np.where(((az1_begin < az_sta_begin)&(az_sta_begin < 360)) | ((0 < az_sta_begin)&(az_sta_begin < az2_begin)))[0]
        else:
            tmp_idx_begin = np.where((az1_begin < az_sta_begin)&(az_sta_begin < az2_begin))[0]
        if 90 <= daz.az < 270:
            tmp_idx_end = np.where(((az1_end < az_sta_end) & (az_sta_end < az2_end)))[0]
        else:
            tmp_idx_end = np.where(((az1_end < az_sta_end) & (az_sta_end < 360)) | ((0 < az_sta_end) & (az_sta_end < az2_end)))[0]
        sta_daz = distaz(self.stalst[:, 0], self.stalst[:, 1], self.cpara.line[0], self.cpara.line[1])
        sta_dis = sind(daz.baz-sta_daz.az) * sta_daz.degreesToKilometers()
        proj_idx = np.where(np.abs(sta_dis) < width)[0]
        tmp_idx = np.intersect1d(tmp_idx_begin, tmp_idx_end)
        final_idx = np.intersect1d(tmp_idx, proj_idx).astype(int)
        if not final_idx.any():
            self.logger.CCPlog.error('No stations within the profile belt with width of {}'.format(width))
            raise ValueError('Satisfied stations not found')
        return final_idx

    def _pierce_project(self, rfsta):
        rfsta['projlat'] = np.zeros_like(rfsta['piercelat'])
        rfsta['projlon'] = np.zeros_like(rfsta['piercelon'])
        for i, dep in enumerate(self.cpara.depth_axis):
            rfsta['projlat'][:, i], rfsta['projlon'][:, i] = geoproject(rfsta['piercelat'][:, i], rfsta['piercelon'][:, i], *self.cpara.line)

    def stack(self):
        """Stack RFs in bins
        """
        if self.cpara.shape == 'circle' or self.cpara.adaptive: 
            field_lat = 'piercelat'
            field_lon = 'piercelon'
        elif self.cpara.shape == 'rect':
            field_lat = 'projlat'
            field_lon = 'projlon'
        else:
            pass
        for i, bin_info in enumerate(self.bin_loca):
            boot_stack = {}
            bin_mu = np.zeros(self.cpara.stack_range.size)
            bin_ci = np.zeros([self.cpara.stack_range.size, 2])
            bin_count = np.zeros(self.cpara.stack_range.size)
            self.logger.CCPlog.info('{}/{} bin from {:.2f} km at lat: {:.3f} lon: {:.3f}'.format(i + 1, self.bin_loca.shape[0], self.profile_range[i], bin_info[0], bin_info[1]))
            if self.cpara.width is None and self.cpara.shape == 'circle':
                idxs = self.idxs[i]
            else:
                idxs = self.idxs
            for j, dep in enumerate(self.cpara.stack_range):
                idx = int(j * self.stack_mul + self.cpara.stack_range[0]/self.cpara.dep_val)
                bin_dep_amp = np.array([])
                for k in idxs:
                    stop_idx = np.where(self.rfdep[k]['stopindex'] >= idx)[0]
                    fall_idx = np.where(distaz(self.rfdep[k][field_lat][stop_idx, idx], self.rfdep[k][field_lon][stop_idx, idx],
                                        bin_info[0], bin_info[1]).delta < self.fzone[j])[0]
                    bin_dep_amp = np.append(bin_dep_amp, self.rfdep[k]['moveout_correct'][stop_idx[fall_idx], idx])
                bin_mu[j], bin_ci[j], bin_count[j] = boot_bin_stack(bin_dep_amp, n_samples=self.cpara.boot_samples)
            boot_stack['bin_lat'] = bin_info[0]
            boot_stack['bin_lon'] = bin_info[1]
            boot_stack['profile_dis'] = self.profile_range[i]
            boot_stack['mu'] = bin_mu
            boot_stack['ci'] = bin_ci
            boot_stack['count'] = bin_count
            self.stack_data.append(boot_stack)   

    def save_stack_data(self, format='npz'):
        """If format is \'npz\', saving stacked data and parameters to local as a npz file. To load the file, please use data = np.load(fname, allow_pickle=True).
        data['cpara'] is the parameters when CCP stacking.
        data['stack_data'] is the result of stacked data.

        If format is \'dat\' the stacked data will be save into a txt file with 8 columns, including bin_lat, bin_lon, profile_dis, depth, amp, ci_low, ci_high and count.
        where bin_lat and bin_lon represent the position of each bin; profile_dis represents the distance in km between each bin and the start point of the profile; depth represents depth of each bin; amp means the stacked amplitude; ci_low and ci_high mean confidence interval with bootstrap method; count represents stacking number of each bin.

        :param format: Format for stacked data
        :type format: str
        """
        self.cpara.stackfile = fix_filename(self.cpara.stackfile, format)
        self.logger.CCPlog.info('Saving stacked data to {}'.format(self.cpara.stackfile))
        if not isinstance(self.cpara.stackfile, str):
            self.logger.CCPlog.error('fname should be in \'str\'')
            raise ValueError('fname should be in \'str\'')
        if format == 'npz':
            np.savez(self.cpara.stackfile, cpara=self.cpara, stack_data=self.stack_data)
        elif format == 'dat':
            with open(self.cpara.stackfile, 'w') as f:
                for i, bin in enumerate(self.stack_data):
                    for j, dep in enumerate(self.cpara.stack_range):
                        if dep == self.cpara.stack_range[-1]:
                            f.write('{:.4f}\t{:.4f}\t{:.4f}\t{:.2f} nan nan nan nan\n'.format(bin['bin_lat'], bin['bin_lon'],
                                                                                              self.profile_range[i], dep))
                        else:
                            f.write('{:.4f}\t{:.4f}\t{:.4f}\t{:.2f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:d}\n'.format(bin['bin_lat'], bin['bin_lon'],
                                                                                    self.profile_range[i],
                                                                                    dep, bin['mu'][j], bin['ci'][j, 0],
                                                                                    bin['ci'][j, 1], int(bin['count'][j])))


if __name__ == '__main__':
    bin_loca, _ = init_profile(27.5, 94, 36.5, 92, 5)
    print(bin_loca.shape)
