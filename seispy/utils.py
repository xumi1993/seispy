from os.path import join, dirname, exists, abspath
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
from seispy import geo
from seispy.geo import geo2sph, km2deg, skm2srad, sph2geo, srad2skm
from seispy import distaz
import numpy as np
from scipy.interpolate import interp1d, interpn
import seispy


def load_cyan_map():
    path = join(dirname(__file__), 'data', 'cyan.mat')
    carray = loadmat(path)['cyan']
    return ListedColormap(carray)


def check_path(key, path):
    if not exists(path):
        raise FileNotFoundError('No such file or directory of {}: {}'.format(key, path))
    else:
        return path

def check_stack_val(stack_val, dep_val):
    if np.mod(stack_val, dep_val) == 0:
        return stack_val/dep_val
    else:
        raise ValueError('stack_val must be a multiple of dep_val')


def read_rfdep(path):
    try:
        return np.load(path, allow_pickle=True)
    except:
        try:
            return np.load(path+'.npy', allow_pickle=True)
        except Exception as e:
            raise FileNotFoundError('Cannot open file of {}'.format(path))


class DepModel(object):
    def __init__(self, YAxisRange, velmod='iasp91', elevation=0):
        self.elevation = elevation
        VelocityModel = np.loadtxt(self.from_file(velmod))
        self.depthsraw = VelocityModel[:, 0]
        self.vpraw = VelocityModel[:, 1]
        self.vsraw = VelocityModel[:, 2]
        self.depths = YAxisRange.astype(float)
        self.dep_val = np.average(np.diff(self.depths))
        if elevation == 0:
            self.depths_elev = self.depths
            self.depths_extend = self.depths
        else:
            dep_append = np.arange(self.depths[-1]+self.dep_val, 
                               self.depths[-1]+self.dep_val+np.floor(elevation/self.dep_val+1), self.dep_val)
            self.depths_extend = np.append(self.depths, dep_append)
            self.depths_elev = np.append(self.depths, dep_append) - elevation
        self.dz = np.append(0, np.diff(self.depths_extend))
        self.vp = interp1d(self.depthsraw, self.vpraw, bounds_error=False,
                           fill_value=self.vpraw[0])(self.depths_elev)
        self.vs = interp1d(self.depthsraw, self.vsraw, bounds_error=False,
                           fill_value=self.vsraw[0])(self.depths_elev)
        self.R = 6371.0 - self.depths_elev

    def from_file(self, mode_name):
        if exists(mode_name):
            filename = mode_name
        elif exists(join(dirname(__file__), 'data', mode_name.lower()+'.vel')):
            filename = join(dirname(__file__), 'data', mode_name.lower()+'.vel')
        else:
            raise ValueError('No such file of velocity model')
        return filename

    def tpds(self, rayps, raypp, sphere=True):
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        tps = np.cumsum((np.sqrt((radius / self.vs) ** 2 - rayps ** 2) -
                        np.sqrt((radius / self.vp) ** 2 - raypp ** 2)) *
                        (self.dz / radius))
        return tps
    
    def tpppds(self, rayps, raypp, sphere=True):
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        tps = np.cumsum((np.sqrt((radius / self.vs) ** 2 - rayps ** 2) +
                        np.sqrt((radius / self.vp) ** 2 - raypp ** 2)) *
                        (self.dz / radius))
        return tps
    
    def tpspds(self, rayps, sphere=True):
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        tps = np.cumsum(2*np.sqrt((radius / self.vs) ** 2 - rayps ** 2)*
                        (self.dz / radius))
        return tps

    def radius_s(self, rayp, phase='P', sphere=True):
        if phase == 'P':
            vel = self.vp
        else:
            vel = self.vs
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        hor_dis = np.cumsum((self.dz / radius) / np.sqrt((1. / (rayp ** 2. * (radius / vel) ** -2)) - 1))
        return hor_dis

    def raylength(self, rayp, phase='P', sphere=True):
        if phase == 'P':
            vel = self.vp
        else:
            vel = self.vs
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        raylen = (self.dz * radius) / (np.sqrt(((radius / self.vs) ** 2) - (rayp ** 2)) * vel)
        return raylen


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


def create_center_bin_profile(stations, val=5, method='linear'):
    if not isinstance(stations, seispy.rf2depth_makedata.Station):
        raise TypeError('Stations should be seispy.rf2depth_makedata.Station')
    dis_sta = prof_range(stations.stla, stations.stlo)
    dis_inter = np.append(np.arange(0, dis_sta[-1], val), dis_sta[-1])
    r, theta, phi = geo2sph(np.zeros(stations.stla.size), stations.stla, stations.stlo)
    # t_po = np.arange(stations.stla.size)
    # ip_t_po = np.linspace(0, stations.stla.size, bin_num)
    theta_i = interp1d(dis_sta, theta, kind=method, bounds_error=False, fill_value='extrapolate')(dis_inter)
    phi_i = interp1d(dis_sta, phi, kind=method, bounds_error=False, fill_value='extrapolate')(dis_inter)
    _, lat, lon = sph2geo(r, theta_i, phi_i)
    # dis = prof_range(lat, lon)
    return lat, lon, dis_inter


def prof_range(lat, lon):
    dis = [0]
    for i in range(lat.size-1):
        dis.append(distaz(lat[i], lon[i], lat[i+1], lon[i+1]).degreesToKilometers())
    return np.cumsum(dis)