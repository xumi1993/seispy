from os.path import join, dirname, exists
from scipy.io import loadmat
from matplotlib.colors import ListedColormap


def load_cyan_map():
    path = join(dirname(__file__), 'data', 'cyan.mat')
    carray = loadmat(path)['cyan']
    return ListedColormap(carray)


def check_path(key, path):
    if not exists(path):
        raise FileNotFoundError('No such file or directory of {}: {}'.format(key, path))
    else:
        return path