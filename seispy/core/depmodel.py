import glob
from os.path import exists, join, dirname
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d, interpn
from seispy.utils import vs2vprho
import warnings

warnings.filterwarnings("ignore", "invalid value encountered in sqrt")


def _search_vel_file(mode_name):
    """
    search vel file given by mode_name
        and try to open it with np.loadtxt:

    1. precise name of vel mod file
    2. default vel mod (containing iasp01
        default vel mod stored in seispy/data/modname.vel file

    do basic value evaluation, cal rho if not given in files
    Parameters
    ----------
    mode_name

    Returns
    ----------
    model matrix [layers:4]:
    depth   vp   vs   rho
    ....    ...  ...  ....
    """
    if not isinstance(mode_name, str):
        raise TypeError('velmod should be in str type')
    if exists(mode_name):
        filename = mode_name
    ## if not found, search in default vel model fold
    elif exists(join(dirname(dirname(__file__)), 'data', mode_name.lower() + '.vel')):
        filename = join(dirname(dirname(__file__)), 'data', mode_name.lower() + '.vel')
    else:
        raise ValueError('No such file of velocity model')

    try:
        raw_model = np.loadtxt(filename)
    except:
        raise IOError
    # check cols of file
    if raw_model.shape[1] < 3 or raw_model.shape[1] > 4:
        raise ValueError('The file should contain 3 or 4 columns')

    # cal rho if rho is not given in vel files.
    if raw_model.shape[1] == 3:
        model = np.zeros((raw_model.shape[0], 4))
        model[:, :3] = raw_model[:, :]
        _p, rho = vs2vprho(raw_model[:, 2])
        model[:, 3] = rho
        return model
    else:
        return raw_model


def _layer2grid(dep_range, model):
    """
    trans model from layer_model to layer_grid

    leave stuffing and interp 2 discretize

    dep_range : grids at depth axis, np.ndarray
    h : thichness of each layer
    vp, vs, rho

    Returns
    -------

    """

    neo_model = np.zeros((len(dep_range), 4)).astype(float)
    picks = np.searchsorted(model[:, 0], dep_range, side="left")
    for _i, _j in enumerate(picks):
        neo_model[_i, :] = model[_j, :]
    return neo_model[:, 1], neo_model[:, 2], neo_model[:, 3]


def _intep_mod(model, depths_elev):
    vp = interp1d(model[:,0], model[:,1], bounds_error=False,
                       fill_value=model[0,1])(depths_elev)
    vs = interp1d(model[:,0], model[:,2], bounds_error=False,
                       fill_value=model[0,2])(depths_elev)
    rho = interp1d(model[:,0], model[:,3], bounds_error=False,
                        fill_value=model[0,3])(depths_elev)
    return vp, vs, rho

def _from_layer_model(dep_range, h, vs, vp=None, rho=None):
    dep = 0
    vp_dep = np.zeros_like(dep_range).astype(float)
    vs_dep = np.zeros_like(dep_range).astype(float)
    rho_dep = np.zeros_like(dep_range).astype(float)
    for i, layer in enumerate(h):
        if (layer == 0 and i == len(h)-1) or \
           (dep+layer < dep_range[-1] and i == len(h)-1):
           idx = np.where(dep_range>=dep)[0]
        else:
            idx = np.where((dep_range >= dep) & (dep_range < dep+layer))[0]
        if idx.size == 0:
            raise ValueError('The thickness of layer {} less than the depth interval'.format(i+1))
        vs_dep[idx] = vs[i]
        if vp is not None:
            vp_dep[idx] = vp[i]
        else:
            vp_dep[idx], _ = vs2vprho(vs[i])
        if rho is not None:
            rho_dep[idx] = rho[i]
        else:
            _, rho_dep[idx] = vs2vprho(vs[i])
        dep += layer
        if dep > dep_range[-1]:
            break
    return np.vstack((dep_range, vp_dep, vs_dep, rho_dep)).T


class DepModel(object):
    """
    Class for constructing 1D velocity model and computing delay time of Ps, PpPs and PsPs Rays.

    .. rubric:: Examples

    >>> model = DepModel(np.array([0, 20.1, 35.1, 100]))
    >>> print(model.dz)
    [ 0.  20.1 15.  64.9]
    >>> print(model.vp)
    [5.8        6.5        8.04001059 8.04764706]
    >>> print(model.vs)
    [3.36       3.75       4.47003177 4.49294118]
    """

    def __init__(self, dep_range, velmod='iasp91', elevation=0., layer_mod=False):
        """ Initialize DepModel object
        
        :type dep_range: numpy.ndarray
        :param dep_range: Depth range for conversion
        :type velmod: str, optional
        :param velmod: Text file of 1D velocity model with first 3 columns of depth/thickness, Vp and Vs, by default 'iasp91'
        :type elevation: float, optional
        :param elevation: Elevation in km, by default 0.0
        :type layer_mod: bool, optional
        :param layer_mod: True for search, and False for interp1d, by default False
        
        """
        self.isrho = False
        self.elevation = elevation
        self.layer_mod = layer_mod
        # dep layer for CCP or other purpose, relative to sea level, dont contain any depth infomation
        self.depths = dep_range.astype(float)
        self.dep_val = np.average(np.diff(self.depths))

        try:
            self.model_array = _search_vel_file(velmod)
        except (IOError, TypeError):
            return
        else:
            self._elevation()
            if layer_mod:
                self.vp, self.vs, self.rho = \
                    _layer2grid(self.depths, self.model_array)
            else:
                self.vp, self.vs, self.rho = \
                    _intep_mod(self.model_array, self.depths_elev)

    @classmethod
    def read_layer_model(cls, dep_range, h, vp, vs, rho=None, elevation=0):
        """ Read layer model from given parameters

        :type dep_range: numpy.ndarray
        :param dep_range: Depth range for conversion
        :type h: numpy.ndarray
        :param h: Thickness of each layer
        :type vp: numpy.ndarray
        :param vp: P-wave velocity of each layer
        :type vs: numpy.ndarray
        :param vs: S-wave velocity of each layer
        :type rho: numpy.ndarray
        :param rho: Density of each layer, by default None
        :type elevation: float
        :param elevation: Elevation in km, by default 0
        :rtype: DepModel
        :return: DepModel object

        .. rubric:: Examples
        
        >>> dep_range = np.arange(100)
        >>> h = np.array([20, 15., 0])
        >>> vp = np.array([5.8, 6.5, 8.04])
        >>> vs = np.array([3.36, 3.75, 4.47])
        >>> model = DepModel.read_layer_model(dep_range, h, vp, vs)
        """
        mod = cls(dep_range, velmod=None, layer_mod=True, elevation=elevation)
        if rho is not None:
            mod.isrho = True
        mod._elevation()
        mod.model_array = _from_layer_model(mod.depths, h, vp, vs, rho=rho)
        mod.vp, mod.vs, mod.rho = _layer2grid(mod.depths, mod.model_array)
        return mod

    def _elevation(self):
        """
        set all depth related values:
        1. depths_elev: depth array contains layer above sea level(represent by value lower than 0
        2. depths_extend: depth array contains layer above sea level( evaluate from 0 to original depth
        3. dz: 0, thick_1, thick_2....
        4. thickness: thick_1, thick_2, ... , 0
        requires elevation, depths_range or other component

        .. rubric:: Examples


        >>> model = DepModel(np.array([0, 20.1, 35.1, 100]))
        >>> model.depths_elev
        array([  0. ,  20.1,  35.1, 100. ])
        >>> model.depths_extend
        array([  0. ,  20.1,  35.1, 100. ])
        >>> model.dz
        array([ 0. , 20.1, 15. , 64.9])
        >>> model.thickness
        array([20.1, 15. , 64.9,  0. ])

        >>> model = DepModel(np.array([0, 20.1, 35.1, 100]),elevation=10.)
        >>> print(model.depths_elev)
        [-10.          10.1         25.1         90.         123.33333333]
        >>> print(model.depths_extend)
        [  0.          20.1         35.1        100.         133.33333333]
        >>> print(model.dz)
        [ 0.         20.1        15.         64.9        33.33333333]
        >>> print(model.thickness)
        [20.1        15.         64.9        33.33333333  0.        ]
        """
        if self.elevation == 0:
            self.depths_elev = self.depths
            self.depths_extend = self.depths
        else:
            depths_append = np.arange(self.depths[-1]+self.dep_val,
                               self.depths[-1]+self.dep_val+np.floor(self.elevation/self.dep_val+1), self.dep_val)

            self.depths_extend = np.append(self.depths, depths_append)
            self.depths_elev = np.append(self.depths, depths_append) - self.elevation
        self.dz = np.append(0, np.diff(self.depths_extend))
        self.thickness = np.append(np.diff(self.depths_extend), 0.)
        self.R = 6371.0 - self.depths_elev

    def plot_model(self, show=True):
        """
        plot model with matplotlib
        
        :type show: bool, optional
        :param show: whether to show the plot, by default True
        """        
        plt.style.use('bmh')
        if self.isrho:
            self.model_fig = plt.figure(figsize=(6, 6))
            fignum = 2
        else:
            self.model_fig = plt.figure(figsize=(4, 6))
            fignum = 1
        self.model_ax = self.model_fig.add_subplot(1, fignum, 1)
        self.model_ax.step(self.vp, self.depths, where='pre', label='Vp')
        self.model_ax.step(self.vs, self.depths, where='pre', label='Vs')
        self.model_ax.legend()
        self.model_ax.set_xlabel('Velocity (km/s)')
        self.model_ax.set_ylabel('Depth (km)')
        self.model_ax.set_ylim([self.depths[0], self.depths[-1]])
        self.model_ax.invert_yaxis()
        if self.isrho:
            self.rho_ax = self.model_fig.add_subplot(1, fignum, 2)
            self.rho_ax.step(self.rho, self.depths, where='pre', color='C2', label='Density')
            self.rho_ax.legend()
            self.rho_ax.set_xlabel('Density (km/s)')
            self.rho_ax.set_ylim([self.depths[0], self.depths[-1]])
            self.rho_ax.invert_yaxis()
        if show:
            plt.show()

    def tpds(self, rayps, raypp, sphere=True):
        # generate docstring
        """
        calculate travel time of Pds

        :type rayps: float or numpy.ndarray
        :param rayps: ray parameter of Ps wave
        :type raypp: float or numpy.ndarray
        :param raypp: ray parameter of P wave
        :type sphere: bool, optional
        :param sphere: whether to use sphere earth, by default True
        :rtype: numpy.ndarray
        :return: travel time of Pds
        """

        if sphere:
            radius = self.R
        else:
            radius = 6371.
        tps = np.cumsum((np.sqrt((radius / self.vs) ** 2 - rayps ** 2) -
                         np.sqrt((radius / self.vp) ** 2 - raypp ** 2)) *
                        (self.dz / radius))
        return tps

    def tpppds(self, rayps, raypp, sphere=True):
        """
        calculate travel time of Ppds

        :type rayps: float or numpy.ndarray
        :param rayps: ray parameter of Ps wave
        :type raypp: float or numpy.ndarray
        :param raypp: ray parameter of P wave
        :type sphere: bool, optional
        :param sphere: whether to use sphere earth, by default True
        :rtype: numpy.ndarray
        :return: travel time of Ppds
        """
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        tps = np.cumsum((np.sqrt((radius / self.vs) ** 2 - rayps ** 2) +
                         np.sqrt((radius / self.vp) ** 2 - raypp ** 2)) *
                        (self.dz / radius))
        return tps

    def tpspds(self, rayps, sphere=True):
        """
        calculate travel time of Ppsds
        
        :type rayps: float or numpy.ndarray
        :param rayps: ray parameter of Ps wave
        :type sphere: bool, optional
        :param sphere: whether to use sphere earth, by default True
        :rtype: numpy.ndarray
        :return: travel time of Ppsds
        """
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        tps = np.cumsum(2 * np.sqrt((radius / self.vs) ** 2 - rayps ** 2) *
                        (self.dz / radius))
        return tps

    def radius_s(self, rayp, phase='P', sphere=True):
        """
        calculate horizontal radius from piercing point to station postion.
        
        P for Sp and S for Ps.

        :type rayp: float or numpy.ndarray
        :param rayp: ray parameter
        :type phase: str, optional
        :param phase: phase name, by default 'P'
        :type sphere: bool, optional
        :param sphere: whether to use sphere earth, by default True
        :rtype: numpy.ndarray
        :return: horizontal radius
        
        .. rubric:: Examples

        >>> model = DepModel(np.array([0, 20.1, 35.1, 100]))
        >>> model.dz
        array([ 0. , 20.1, 15. , 64.9])
        >>> model.R
        array([6371. , 6350.9, 6335.9, 6271. ])
        >>> model.radius_s(1.2,phase="S", sphere=False)*111.2
        array([0.        , 0.0002478 , 0.00046823, 0.00142685])
        """
        if phase == 'P':
            vel = self.vp
        else:
            vel = self.vs
        if sphere:
            radius = self.R
        else:
            radius = 6371.
        hor_dis = np.cumsum((self.dz / radius) / np.sqrt((1. / (rayp ** 2. * (radius / vel) ** -2)) - 1))
        #hor_dis = np.sqrt((1. / (rayp ** 2. * (radius / vel) ** -2)) - 1)
        return hor_dis

    def save_tvel(self, filename):
        """
        save vel mod in tvel format for taup

        :type filename: str
        :param filename: output file name

        .. rubric:: Examples

        >>> model = DepModel(np.array([0, 20.1, 35.1, 100]))
        >>> model.save_tvel("test")
            0.00     5.80     3.36     2.72
           20.10    5.800     3.36     2.72
           20.10     6.50     3.75     2.92
           35.10    6.500     3.75     2.92
           35.10     8.04     4.47     3.32
          100.00    8.040     4.47     3.32
          100.00     8.05     4.49     3.36
        """
        outlines =[]
        depth = self.depths
        if depth[0] != 0:
            depth[0] = 0
        for _i in range(self.depths.shape[0]):
            out = f'{depth[_i]:>8.2f} {self.vp[_i]:>8.2f} {self.vs[_i]:>8.2f} {self.rho[_i]:>8.2f}'
            outlines.append(out)
            if _i == self.depths.shape[0] - 1:
                break
            out = f'{depth[_i+1]:>8.2f} {self.vp[_i]:>8.3f} {self.vs[_i]:>8.2f} {self.rho[_i]:>8.2f}'
            outlines.append(out)
        # this switch for doctest
        if filename == 'test':
            for _i in range(len(outlines)):
                print(outlines[_i])
        else:
            try:
                f = open(filename, 'w')
                f.write("temp_mod-P\ntemp_mod-S")
                for i in range(len(outlines)):
                    f.write(outlines[i])
                f.close()
            except IOError:
                raise IOError('cannot write to {}'.format(filename))


    def raylength(self, rayp, phase='P', sphere=True):
        """
        calculate ray length, P for Sp and S for Ps
        
        :type rayp: float or numpy.ndarray
        :param rayp: ray parameter
        :type phase: str, optional
        :param phase: phase name, by default 'P'
        :type sphere: bool, optional
        :param sphere: whether to use sphere earth, by default True
        :rtype: numpy.ndarray
        :return: ray length
        """

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

    @classmethod
    def ccp_model(cls, dep_range = np.array([0, 20.1, 35.1, 100]),
                  elevation = 0, layerd = False, **kwargs):
        """
        import ccp configure and init DepModel object for time2depth convertion
        if any parameters given is wrong, return a default DepModel object

        there's 3 types of input is allowed:

        1. mod3d, stla, stlo: for 3d model( need modification

        2. modfolder, staname: for dir

        3. mod: for single file
        """
        if kwargs.get("mod3d", None):
            stla = kwargs.pop("stla",None)
            stlo = kwargs.pop("stlo", None)
            # if stla and stlo is not given return iasp91 model
            if not stla or not stlo:
                return cls(dep_range, "iasp91", elevation, layerd)
            try:
                mod = np.load(kwargs["mod3d"])
                dep_tmp = DepModel(dep_range,"iasp91", elevation, layerd)
                dep_tmp.vp, dep_tmp.vs = interp_depth_model(mod, stla, stlo, dep_range)
                return dep_tmp
            except:
                return cls(dep_range, "iasp91", elevation, layerd)

        elif kwargs.get("modfolder", None):
            try:
                mod = _load_mod(kwargs.pop("modfolder", "./"), kwargs.pop("staname", None))
            except:
                mod = "iasp91"
            finally:
                return cls(dep_range, mod, elevation, layerd)
        else:
            mod = kwargs.pop("mod", "iasp91")
            return cls(dep_range,mod, elevation, layerd)

def _load_mod(datapath, staname):
    """Load 1D velocity model files with suffix of ".vel". The model file should be including 3 columns with depth, vp and vs.

    :param datapath: Folder name with 1D velocity model files.
    :type datapath: string
    :param staname: The station name as a part of file name of 1D velocity model files.
    :type staname: string
    """
    expresion = join(datapath, "*"+staname+"*.vel")
    modfiles = glob.glob(expresion)
    if len(modfiles) == 0:
        raise FileNotFoundError("The model file of {} were not found.".format(expresion))
    elif len(modfiles) > 1:
        raise ValueError('More then 1 file were found as the expresion: {}'.format(expresion))
    else:
        return modfiles[0]

def interp_depth_model(model, lat, lon, new_dep):
    """ Interpolate Vp and Vs from 3D velocity with a specified depth range.

    :param mod3d: 3D velocity loaded from a ``.npz`` file
    :type mod3d: :meth:`np.lib.npyio.NpzFile`
    :param lat: Latitude of position in 3D velocity model
    :type lat: float
    :param lon: Longitude of position in 3D velocity model
    :type lon: float
    :param new_dep: 1D array of depths in km
    :type new_dep: :meth:`np.ndarray`
    :rtype: :meth:`np.ndarray`
    :return: Vp and Vs in ``new_dep``
    
    """
    #  model = np.load(modpath)
    points = [[depth, lat, lon] for depth in new_dep]
    vp = interpn((model['dep'], model['lat'], model['lon']), model['vp'], points, bounds_error=False, fill_value=None)
    vs = interpn((model['dep'], model['lat'], model['lon']), model['vs'], points, bounds_error=False, fill_value=None)
    return vp, vs

if __name__ == "__main__":
    import doctest
    doctest.testmod()
