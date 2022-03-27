import numpy as np
from seispy.geo import km2deg, latlon_from, cosd, extrema, skm2srad, rad2deg
from seispy import distaz
from seispy.rfcorrect import DepModel
from seispy.setuplog import setuplog
from scikits.bootstrap import ci
from seispy.ccppara import ccppara, CCPPara
from seispy.signal import smooth
from seispy.utils import check_stack_val, read_rfdep
import warnings


def gen_center_bin(center_lat, center_lon, len_lat, len_lon, val):
    """
Create spaced grid point with coordinates of the center point in the area in spherical coordinates.

:: 
    ---------------------------------------------------------
    |                           |                           |
    |                           |                           |
    |                        len_lon                        |
    |                           |                           |
    |                           |                           |
    ---- len_lat --- (center_lon, center_lat) --- len_lat ---
    |                           |                           |
    |                           |                           |
    |                        len_lon                        |
    |                           |                           |
    |                           |                           |
    ---------------------------------------------------------

:param center_lat: Latitude of the center point.
:type center_lat: float
:param center_lon: Longitude of the center point.
:type center_lon: float
:param len_lat: Half length in degree along latitude axis.
:type len_lat: float
:param len_lon: Half length in degree along longitude axis.
:type len_lon: float
:param val: Interval in degree between adjacent grid point.
:type val: float
:return: Coordinates of Grid points.
:rtype: 2-D ndarray of floats with shape (n, 2), where n is the number of grid points.
    """
    lats = np.arange(0, 2*len_lat, val)
    lons = np.arange(0, 2*len_lon, val)
    plat, plon = latlon_from(center_lat, center_lon, 0, 90)
    da = distaz(plat, plon, center_lat, center_lon)
    begx = -len_lon 
    begy = -len_lat
    bin_loca = []
    for j in range(lats.size):
        delyinc = j * val + begy
        delt = da.delta + delyinc
        for i in range(lons.size):
            azim = da.az + (begx + i * val) / cosd(delyinc)
            glat, glon = latlon_from(plat, plon, azim, delt)
            if glon > 180:
                glon -= 360
            bin_loca.append([glat, glon])
    return np.array(bin_loca)


def bin_shape(cpara):
    if cpara.shape == 'rect':
        raise ValueError('The shape of bins must be set to \'circle\' in ccp3d mode.')
    if cpara.bin_radius is None:
        depmod = DepModel(cpara.stack_range)
        fzone = km2deg(np.sqrt(0.5*cpara.domperiod*depmod.vs*cpara.stack_range))
    else:
        fzone = np.ones_like(cpara.stack_range) * km2deg(cpara.bin_radius)
    return fzone


def boot_bin_stack(data_bin, n_samples=3000):
    warnings.filterwarnings("ignore")
    count = data_bin.shape[0]
    if count > 1:
        if n_samples is not None:
            cci = ci(data_bin, n_samples=n_samples)
        else:
            cci = np.array([np.nan, np.nan])
        mu = np.average(data_bin)
    else:
        cci = np.array([np.nan, np.nan])
        mu = np.nan
    return mu, cci, count


def _get_sta(rfdep):
    return np.array([[sta['stalat'], sta['stalon']] for sta in rfdep])


def _sta_val(stack_range, radius):
    dep_mod = DepModel(stack_range)
    x_s = np.cumsum((dep_mod.dz / dep_mod.R) / np.sqrt((1. / (skm2srad(0.08) ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))
    dis = radius + rad2deg(x_s[-1]) + 0.5
    return dis


class CCP3D():
    def __init__(self, cfg_file=None, log=None):
        """Class for 3-D CCP stacking, Usually used to study mantle transition zone structure.

        :param cfg_file: Path to configure file. If not defined a instance of CCP3D.cpara will be initialed, defaults to None
        :type cfg_file: str, optional
        :param log: A logger instance. If not defined, seispy.sutuplog.logger will be initialed, defaults to None
        :type log: seispy.sutuplog.logger , optional
        """
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
            self.logger.CCPlog('Cannot open configure file {}'.format(cfg_file))
            raise FileNotFoundError('{}'.format(e))
        try:
            self.stack_mul = check_stack_val(self.cpara.stack_val, self.cpara.dep_val)
        except Exception as e:
            self.logger.CCPlog.error('{}'.format(e))
            raise ValueError('{}'.format(e))
    
    def read_rfdep(self):
        self.logger.CCPlog.info('Loading RFdepth data from {}'.format(self.cpara.depthdat))
        try:
            self.rfdep = read_rfdep(self.cpara.depthdat)
        except FileNotFoundError as e:
            self.logger.CCPlog.error('{}'.format(e))
            raise FileNotFoundError('Cannot open file of {}'.format(self.cpara.depthdat))

    def initial_grid(self):
        self.read_rfdep()
        self.bin_loca = gen_center_bin(*self.cpara.center_bin)
        self.fzone = bin_shape(self.cpara)
        self.stalst = _get_sta(self.rfdep)
        self.dismin = _sta_val(self.cpara.stack_range, self.fzone[-1])

    def _select_sta(self, bin_lat, bin_lon):
        return np.where(distaz(bin_lat, bin_lon, self.stalst[:, 0], self.stalst[:, 1]).delta <= self.dismin)[0]

    def stack(self):
        """Search conversion points falling within a bin and stack them with bootstrap method.
        """
        for i, bin_info in enumerate(self.bin_loca):
            boot_stack = {}
            bin_mu = np.zeros(self.cpara.stack_range.size)
            bin_ci = np.zeros([self.cpara.stack_range.size, 2])
            bin_count = np.zeros(self.cpara.stack_range.size)
            self.logger.CCPlog.info('{}/{} bin at lat: {:.3f} lon: {:.3f}'.format(i + 1, self.bin_loca.shape[0],
                                                                                  bin_info[0], bin_info[1]))
            idxs = self._select_sta(bin_info[0], bin_info[1])
            for j, dep in enumerate(self.cpara.stack_range):
                idx = int(j * self.stack_mul + self.cpara.stack_range[0]/self.cpara.dep_val)
                bin_dep_amp = np.array([])
                for k in idxs:
                    stop_idx = np.where(self.rfdep[k]['stopindex'] >= idx)[0]
                    fall_idx = np.where(distaz(self.rfdep[k]['piercelat'][stop_idx, idx], self.rfdep[k]['piercelon'][stop_idx, idx],
                                        bin_info[0], bin_info[1]).delta < self.fzone[j])[0]
                    bin_dep_amp = np.append(bin_dep_amp, self.rfdep[k]['moveout_correct'][stop_idx[fall_idx], idx])
                bin_mu[j], bin_ci[j], bin_count[j] = boot_bin_stack(bin_dep_amp, n_samples=self.cpara.boot_samples)
            boot_stack['bin_lat'] = bin_info[0]
            boot_stack['bin_lon'] = bin_info[1]
            boot_stack['mu'] = bin_mu
            boot_stack['ci'] = bin_ci
            boot_stack['count'] = bin_count
            self.stack_data.append(boot_stack)
 
    def save_stack_data(self, fname):
        """Save stacked data and parameters to local as a npz file. To load the file, please use data = np.load(fname, allow_pickle=True).
        data['cpara'] is the parameters when CCP stacking.
        data['stack_data'] is the result of stacked data.

        :param fname: file name of stacked data
        :type fname: str
        """
        if not isinstance(fname, str):
            self.logger.CCPlog.error('fname should be in \'str\'')
            raise ValueError('fname should be in \'str\'')
        np.savez(fname, cpara=self.cpara, stack_data=self.stack_data)
    
    def _search_peak(self, tr, peak_410_min=380, peak_410_max=440, peak_660_min=630, peak_660_max=690):
        tr = smooth(tr, half_len=4)
        idx_all = extrema(tr)
        idx_ex = np.where(tr[idx_all] > 0)[0]
        #idx = idx_all[idx_ex]
        peak = tr[idx_all[idx_ex]]
        peak_depth = idx_all[idx_ex] * self.cpara.stack_val + self.cpara.stack_range[0]
        idx_410 = np.where((peak_depth>peak_410_min) & (peak_depth<peak_410_max))[0]
        try:
            idx_410_max = idx_410[np.nanargmax(peak[idx_410])]
            dep_410 = peak_depth[idx_410_max]
        except:
            dep_410 = np.nan
        idx_660 = np.where((peak_depth>peak_660_min) & (peak_depth<peak_660_max))[0]
        try:
            idx_660_max = idx_660[np.nanargmax(peak[idx_660])]
            dep_660 = peak_depth[idx_660_max]
        except:
            dep_660 = np.nan
        return dep_410, dep_660

    def search_good_410_660(self, peak_410_min=380, peak_410_max=440, peak_660_min=630, peak_660_max=690):
        self.good_410_660 = np.zeros_like(self.bin_loca)
        for i, boot_stack in enumerate(self.stack_data):
            self.good_410_660[i, 0], self.good_410_660[i, 1] = self._search_peak(boot_stack['mu'], peak_410_min, peak_410_max, peak_660_min, peak_660_max)

    def save_good_410_660(self, fname):
        with open(fname, 'w') as f:
            for i, good_peak in enumerate(self.good_410_660):
                if np.isnan(good_peak[0]):
                    ci_410 = np.array([np.nan, np.nan])
                    count_410 = np.nan
                else:
                    idx = int((good_peak[0]-self.cpara.stack_range[0]) / self.cpara.stack_val)
                    ci_410 = self.stack_data[i]['ci'][idx]
                    count_410 = self.stack_data[i]['count'][idx]
                if np.isnan(good_peak[1]):
                    ci_660 = np.array([np.nan, np.nan])
                    count_660 = np.nan
                else:
                    idx = int((good_peak[1]-self.cpara.stack_range[0]) / self.cpara.stack_val)
                    ci_660 = self.stack_data[i]['ci'][idx]
                    count_660 = self.stack_data[i]['count'][idx]
                f.write('{:.3f} {:.3f} {:.0f} {:.4f} {:.4f} {:.0f} {:.0f} {:.4f} {:.4f} {:.0f}\n'.format(   
                        self.bin_loca[i, 0], self.bin_loca[i, 1], good_peak[0], ci_410[0], ci_410[1], count_410,
                        good_peak[1], ci_660[0], ci_660[1], count_660))

    @classmethod
    def read_stack_data(cls, stack_data_path, cfg_file=None):
        ccp = cls(cfg_file)
        data = np.load(stack_data_path, allow_pickle=True)
        ccp.stack_data = data['stack_data']
        ccp.cpara = data['cpara'].any()
        ccp.bin_loca = np.array([[bin_info['bin_lat'], bin_info['bin_lon']] for bin_info in ccp.stack_data])
        return ccp


if __name__ == '__main__':
    bin_loca = gen_center_bin(48.5, 100, 5, 8, km2deg(55))
    with open('/workspace/WMHG_MTZ/ccp_results/bin_loca.dat', 'w') as f:
        for binin in bin_loca:
            f.write('{} {}\n'.format(binin[1], binin[0]))
