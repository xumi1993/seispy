import numpy as np
from os.path import expanduser, join, dirname, exists
import configparser
from seispy.geo import km2deg
import warnings


class CCPPara(object):
    def __init__(self):
        self.rfpath = expanduser('~')
        self.rayp_lib = None
        self.depthdat = 'RFdepth.npy'
        self.stackfile = 'ccp.dat'
        self.stalist = 'sta.lst'
        self.peakfile = 'good_410_660.dat'
        self.adaptive= False
        self.velmod = ''
        self.stack_sta_list = ''
        self.domperiod = 5
        self.shape = 'circle'
        self.slide_val = 5
        self.width = 100
        self.bin_radius = 50
        self.line = np.array([])
        self.depth_axis = np.array([])
        self.stack_range = np.array([])
        self.center_bin = []
        self.dep_val = 1
        self.stack_val = 1
        self.boot_samples = None
        self.phase = 1
    
    def __str__(self):
        head = ['{}: {}'.format(k, v) for k, v in self.__dict__.items()]
        return '\n'.join(head)

    @property
    def bin_radius(self):
        return self._bin_radius

    @bin_radius.setter
    def bin_radius(self, value):
        if not (isinstance(value, (int, float)) or value is None):
            raise TypeError('Error type of bin_radius')
        else:
            self._bin_radius = value

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, value):
        if not isinstance(value, str):
            raise TypeError('ccppara.shape must be str type')
        elif value.lower() not in ('circle', 'rect'):
            raise ValueError('ccppara.shape must be in \'circle\' or \'rect\'')
        else:
            self._shape = value.lower()


def ccppara(cfg_file):
    cpara = CCPPara()
    cf = configparser.ConfigParser()
    try:
        cf.read(cfg_file)
    except Exception:
        raise FileNotFoundError('Cannot open configure file %s' % cfg_file)
    # para for FileIO section
    cpara.rfpath = cf.get('FileIO', 'rfpath')
    rayp_lib = cf.get('FileIO', 'rayp_lib')
    if rayp_lib == '':
        cpara.rayp_lib = None
    else:
        cpara.rayp_lib = rayp_lib
    cpara.depthdat = cf.get('FileIO', 'depthdat')
    cpara.stackfile = cf.get('FileIO', 'stackfile')
    cpara.stalist = cf.get('FileIO', 'stalist')
    cpara.stack_sta_list = cf.get('FileIO', 'stack_sta_list')
    if cf.has_option('FileIO', 'peakfile'):
        fname = cf.get('FileIO', 'peakfile')
        if fname != '':
            cpara.peakfile = fname
    velmod = cf.get('FileIO', 'velmod')
    if velmod == '':
        cpara.velmod = join(dirname(__file__), 'data', 'iasp91.vel')
    elif not exists(velmod):
        cpara.velmod = join(dirname(__file__), 'data', '{}.vel'.format(velmod.lower()))
    else:
        cpara.velmod = velmod
    # para for bin section
    cpara.shape = cf.get('bin', 'shape')
    try:
        cpara.adaptive = cf.getboolean('bin', 'adaptive')
    except:
        cpara.adaptive = False
    try:
        cpara.domperiod = cf.getfloat('bin', 'domperiod')
    except:
        cpara.domperiod = None
    try:
        cpara.width = cf.getfloat('bin', 'width')
    except:
        cpara.width = None
    try:
        cpara.bin_radius = cf.getfloat('bin', 'bin_radius')
    except:
        cpara.bin_radius = None
    try:
        cpara.slide_val = cf.getfloat('bin', 'slide_val')
    except:
        warnings.warn('slide_val not found. Setup it for CCP stacking')
    try:
        # para for line section
        lat1 = cf.getfloat('line', 'profile_lat1')
        lon1 = cf.getfloat('line', 'profile_lon1')
        lat2 = cf.getfloat('line', 'profile_lat2')
        lon2 = cf.getfloat('line', 'profile_lon2')
        cpara.line = np.array([lat1, lon1, lat2, lon2])
    except:
        warnings.warn('line section not found. Setup it for ccp_profile')
    try:
        # para for center bins
        cla = cf.getfloat('spacedbins', 'center_lat')
        clo = cf.getfloat('spacedbins', 'center_lon')
        hlla = cf.getfloat('spacedbins', 'half_len_lat')
        hllo = cf.getfloat('spacedbins', 'half_len_lon')
        cpara.center_bin = [cla, clo, hlla, hllo, km2deg(cpara.slide_val)]
    except:
        warnings.warn('No such section of spacedbins and line. Setup them for CCP stacking')
    # para for depth section
    dep_end = cf.getfloat('depth', 'dep_end')
    cpara.dep_val = cf.getfloat('depth', 'dep_val')
    try:
        cpara.phase = cf.getint('depth', 'phase')
    except:
        cpara.phase = 1
    cpara.depth_axis = np.append(np.arange(0, dep_end, cpara.dep_val), dep_end)

    stack_start = cf.getfloat('stack', 'stack_start')
    stack_end = cf.getfloat('stack', 'stack_end')
    cpara.stack_val = cf.getfloat('stack', 'stack_val')
    cpara.stack_range = np.arange(stack_start, stack_end, cpara.stack_val)
    try:
        cpara.boot_samples = cf.getint('stack', 'boot_samples')
    except:
        cpara.boot_samples = None

    return cpara
