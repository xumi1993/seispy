from os.path import join, dirname, exists
import array
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
import numpy as np
from scipy.interpolate import interp1d, interpn
import pandas as pd


def vs2vprho(vs):
    """
    Calculate vp and rho from vs using Brocher (2005) formula.

    :param vs: vs in km/s
    :type vs: float or np.ndarray

    :return: vp and rho in km/s and g/cm^3
    :rtype: tuple
    """
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