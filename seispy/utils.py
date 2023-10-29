from os.path import join, dirname, exists
import array
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
import numpy as np
from scipy.interpolate import interp1d, interpn
import pandas as pd


def vs2vprho(vs):
    vp = 0.9409 + 2.0947*vs - 0.8206*vs**2 + 0.2683*vs**3 - 0.0251*vs**4
    rho = 1.6612*vp - 0.4721*vp**2 + 0.0671*vp**3 - 0.0043*vp**4 + 0.000106*vp**5
    # vs = 0.7858 - 1.2344*vp + 0.7949*vp**2 - 0.1238*vp**3 + 0.0064*vp**4
    return vp, rho


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


def scalar_instance(v):
    if isinstance(
        v,
        (int,
         float,
         np.floating,
         np.integer
        )
    ):
        return True
    else:
        return False


def array_instance(v):
    if isinstance(
        v,
        (list,
         tuple,
         np.ndarray,
         pd.Series,
         array.array,
        )
    ):
        return True
    else:
        return False