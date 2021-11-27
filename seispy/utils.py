from os.path import join, dirname, exists
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
import numpy as np
from scipy.interpolate import interp1d, interpn


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
    def __init__(self, YAxisRange, velmod='iasp91'):
        VelocityModel = np.loadtxt(from_file(velmod))
        self.depthsraw = VelocityModel[:, 0]
        self.vpraw = VelocityModel[:, 1]
        self.vsraw = VelocityModel[:, 2]
        self.vp = interp1d(self.depthsraw, self.vpraw, bounds_error=False, fill_value='extrapolate')(YAxisRange)
        self.vs = interp1d(self.depthsraw, self.vsraw, bounds_error=False, fill_value='extrapolate')(YAxisRange)
        self.depths = YAxisRange
        self.dz = np.append(0, np.diff(self.depths))
        self.R = 6371 - self.depths


def from_file(mode_name):
    if exists(mode_name):
        filename = mode_name
    elif exists(join(dirname(__file__), 'data', mode_name.lower()+'.vel')):
        filename = join(dirname(__file__), 'data', mode_name.lower()+'.vel')
    else:
        raise ValueError('No such velocity mode')
    return filename


def tpds(dep_mod, rayps, raypp):
    return np.cumsum((np.sqrt((dep_mod.R / dep_mod.vs) ** 2 - rayps ** 2) -
                     np.sqrt((dep_mod.R / dep_mod.vp) ** 2 - raypp ** 2)) *
                     (dep_mod.dz / dep_mod.R))


def radius_s(dep_mod, rayps):
    return np.cumsum((dep_mod.dz / dep_mod.R) /
                     np.sqrt((1. / (rayps ** 2. * (dep_mod.R / dep_mod.vs) ** -2)) - 1))


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