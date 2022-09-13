from os.path import join, dirname, exists, abspath
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
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



def _from_layer_model(h, vp, vs, dep_range):
    dep = 0
    vp_dep = np.zeros_like(dep_range).astype(float)
    vs_dep = np.zeros_like(dep_range).astype(float)
    for i, layer in enumerate(h):
        if (layer == 0 and i == h.size-1) or \
           (dep+layer < dep_range[-1] and i == h.size-1):
           idx = np.where(dep_range>=dep)[0]
        else:
            idx = np.where((dep_range >= dep) & (dep_range < dep+layer))[0]
        if idx.size == 0:
            raise ValueError('The thickness of layer {} less than the depth interval'.format(i+1))
        vp_dep[idx] = vp[i]
        vs_dep[idx] = vs[i]
        dep += layer
        if dep > dep_range[-1]:
            break
    return vp_dep, vs_dep


class DepModel(object):
    def __init__(self, dep_range, velmod='iasp91', elevation=0., layer_mod=False):
        """Class for computing back projection of Ps Ray paths.

        Parameters
        ----------
        dep_range : numpy.ndarray
            Depth range for conversion
        velmod : str, optional
            Text file of 1D velocity model with first 3 columns of depth/thickness, Vp and Vs,
            by default 'iasp91'
        elevation : float, optional
            Elevation in km, by default 0.0
        """
        self.elevation = elevation
        self.layer_mod = layer_mod
        self.depths = dep_range.astype(float)
        try:
            self.model_array = np.loadtxt(self.from_file(velmod))
        except (ValueError, TypeError):
            return
        else:
            self.read_model_file()
            self.discretize()

    def read_model_file(self):
        self.dep_val = np.average(np.diff(self.depths))
        if self.layer_mod:
            self.depthsraw = self.depths
            self.vpraw, self.vsraw = _from_layer_model(self.model_array[:, 0],
                                                       self.model_array[:, 1],
                                                       self.model_array[:, 2],
                                                       self.depths)
        else:
            self.depthsraw = self.model_array[:, 0]
            self.vpraw = self.model_array[:, 1]
            self.vsraw = self.model_array[:, 2]
    
    @classmethod
    def read_layer_model(cls, dep_range, h, vp, vs, **kwargs):
        mod = cls(dep_range, velmod=None, layer_mod=True, **kwargs)
        mod.depthsraw = mod.depths
        mod.vpraw, mod.vsraw = _from_layer_model(h, vp, vs, mod.depths)
        mod.discretize()
        return mod
    
    def plot_model(self, show=True):
        plt.style.use('bmh')
        self.model_fig = plt.figure(figsize=(4,6))
        self.model_ax = self.model_fig.add_subplot()
        self.model_ax.step(self.vp, self.depths, where='post', label='Vp')
        self.model_ax.step(self.vs, self.depths, where='post', label='Vs')
        self.model_ax.legend()
        self.model_ax.set_xlabel('Velocity (km/s)')
        self.model_ax.set_ylabel('Depth (km)')
        self.model_ax.set_ylim([self.depths[0], self.depths[-1]])
        self.model_ax.invert_yaxis()
        if show:
            plt.show()

    def discretize(self):
        if self.elevation == 0:
            self.depths_elev = self.depths
            self.depths_extend = self.depths
        else:
            dep_append = np.arange(self.depths[-1]+self.dep_val, 
                               self.depths[-1]+self.dep_val+np.floor(self.elevation/self.dep_val+1), self.dep_val)
            self.depths_extend = np.append(self.depths, dep_append)
            self.depths_elev = np.append(self.depths, dep_append) - self.elevation
        self.dz = np.append(0, np.diff(self.depths_extend))
        self.vp = interp1d(self.depthsraw, self.vpraw, bounds_error=False,
                           fill_value=self.vpraw[0])(self.depths_elev)
        self.vs = interp1d(self.depthsraw, self.vsraw, bounds_error=False,
                           fill_value=self.vsraw[0])(self.depths_elev)
        self.R = 6371.0 - self.depths_elev

    def from_file(self, mode_name):
        if not isinstance(mode_name, str):
            raise TypeError('velmod should be in str type')
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