from os.path import expanduser
import configparser
import numpy as np


class HKPara(object):
    def __init__(self):
        self.rfpath = expanduser('~')
        self.hkpath = expanduser('~')
        self.hklist = 'hk.dat'
        self.hrange = np.arange(20, 80, 0.1)
        self.krange = np.arange(1.6, 1.9, 0.01)
        self.vp = 6.3

    @property
    def hrange(self):
        return self._hrange

    @hrange.setter
    def hrange(self, value):
        if not (isinstance(value, np.ndarray) or value is None):
            raise TypeError('Error type of hrange')
        else:
            self._hrange = value

    @property
    def krange(self):
        return self._krange

    @krange.setter
    def krange(self, value):
        if not (isinstance(value, np.ndarray) or value is None):
            raise TypeError('Error type of krange')
        else:
            self._krange = value


def hkpara(cfg_file):
    hpara = HKPara()
    cf = configparser.ConfigParser()
    try:
        cf.read(cfg_file)
    except Exception:
        raise FileNotFoundError('Cannot open configure file %s' % cfg_file)

    # para for FileIO section
    hpara.rfpath = cf.get('FileIO', 'rfpath')
    hpara.hkpath = cf.get('FileIO', 'hkpath')
    hpara.hklst = cf.get('FileIO', 'hklst')

    hmin = cf.getfloat('hk', 'hmin')
    hmax = cf.getfloat('hk', 'hmax')
    kmin = cf.getfloat('hk', 'kmin')
    kmax = cf.getfloat('hk', 'kmax')
    hpara.hrange = np.arange(hmin, hmax, 0.1)
    hpara.krange = np.arange(kmin, kmax, 0.01)

    vp = cf.get('hk', 'vp')
    if vp != '':
        hpara.vp = float(vp)
    return hpara
