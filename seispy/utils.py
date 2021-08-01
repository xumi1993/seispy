from os.path import join, dirname
from scipy.io import loadmat
from matplotlib.colors import ListedColormap


def load_cyan_map():
    path = join(dirname(__file__), 'data', 'cyan.mat')
    carray = loadmat(path)['cyan']
    return ListedColormap(carray)
