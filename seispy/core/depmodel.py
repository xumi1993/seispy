from os.path import exists, join, dirname

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from seispy.utils import vs2vprho


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
    elif exists(join(dirname(__file__), '../data', mode_name.lower() + '.vel')):
        filename = join(dirname(__file__), '../data', mode_name.lower() + '.vel')
    else:
        raise ValueError('No such file of velocity model')

    try:
        raw_model = np.loadtxt(filename)
    except:
        raise IOError("failed while reading velocity model {}".format(filename))
    # check cols of file
    if raw_model.shape[1]<3 or raw_model.shape[1]>4:
        raise ValueError("failed while reading velocity model {}, plz check your format".format(filename))

    # cal rho if rho is not given in vel files.
    if raw_model.shape[1]==3:
        model = np.zeros((raw_model.shape[0],4))
        model[:3,:] = raw_model[:,:]
        _p, rho = vs2vprho(raw_model[:,2])
        model[3,:] = rho

        return model
    else:
        return raw_model

def _layer2grid(dep_range, model):
    """
    trans model from layer_model 2 layer_grid
    grid in between discontinuities are set to 0

    leave stuffing and interp 2 discretize

    dep_range : grids at depth axis, np.ndarray
    h : thichness of each layer
    vp, vs, rho

    Returns
    -------

    """

    neo_model = np.zeros((len(dep_range),4)).astype(float)
   # vp_dep = np.zeros_like(dep_range).astype(float)
   # vs_dep = np.zeros_like(dep_range).astype(float)
   # rho_dep = np.zeros_like(dep_range).astype(float)

    picks = np.searchsorted(model[:,0],dep_range, side="left")
    for _i,_j in enumerate(picks):
        neo_model[_i,:]=model[_j,:]
    # return vp_dep, vs_dep, rho_dep
    print(neo_model[:,1])
    return neo_model[:,1], neo_model[:,2], neo_model[:,3]


def _descretize(dep_range, model):
    """


    Parameters
    ----------
    dep_range
    model

    Returns
    -------

    """



def discretize(self,raw_depth, raw_vp, raw_vs, raw_rho):
      """
      full init Depth model by
      interp 1d grid value


      """
      if self.elevation == 0:
          self.depths_elev = self.depths_layer
          self.depths_extend = self.depths_layer
      else:
          dep_append = np.arange(self.depths_layer[-1] + self.dep_val,
                                 self.depths_layer[-1] + self.dep_val + np.floor(self.elevation / self.dep_val + 1), self.dep_val)
          self.depths_extend = np.append(self.depths_layer, dep_append)
          self.depths_elev = np.append(self.depths_layer, dep_append) - self.elevation

      self.dz = np.append(0, np.diff(self.depths_extend))
      self.thickness = np.append(np.diff(self.depths_extend), 0.)

      self.vp = interp1d(raw_depth, raw_vp, bounds_error=False,
                         fill_value=raw_vp[0])(self.depths_elev)
      self.vs = interp1d(raw_depth, raw_vs, bounds_error=False,
                         fill_value=raw_vs[0])(self.depths_elev)
      self.rho = interp1d(raw_depth, raw_vs, bounds_error=False,
                         fill_value=raw_rho[0])(self.depths_elev)

      self.R = 6371.0 - self.depths_elev


def _diff_elev(elevatioon=0, dep_range):
    """
    elevatioon
    dep_range
    ------------------
    gens:
    depths_elev
    depths_extend
    dz
    thickness
    vp, vs, rho,R


    """
    if elevatioon == 0:
        depths_elev = dep_range
        depths_extend = dep_range
    else:
        dep_append = np.arange(self.depths[-1] + self.dep_val, self.elevation + self.dep_val, self.dep_val)
        self.depths_elev = np.concatenate((self.depths - self.elevation, dep_append))

    self.depths_extend = np.concatenate((self.depths, dep_append))
    self.dz = np.diff(self.depths_extend, prepend=0)
    self.thickness = np.diff(self.depths_extend, append=0)


    depths_extend = np.concatenate((depths, dep_append))
    dz = np.diff(depths_extend, prepend=0)
    thickness = np.diff(depths_extend, append=0)



class DepModel(object):
    """
    radiu_s is used to call piercing point
    tps,tpsps or so are used to cal displacement

    """
    def __init__(self, dep_range, velmod='iasp91', elevation=0., layer_mod=False):
        """Class for computing back projection of Ps Ray paths.

        Parameters
        ----------
        dep_range : numpy.ndarray
            Depth range for conversion
            Depth for each layer, not thickness
        velmod : str, optional
            Text file of 1D velocity model with first 3 columns of depth/thickness, Vp and Vs,
            by default 'iasp91'
        elevation : float, optional Elevation in km, by default 0.0
        layer_mod: True for search, and False for interp1d
        """
        self.isrho = False
        self.elevation = elevation
        self.layer_mod = layer_mod
        self.depths_layer = dep_range.astype(float)
        self.dep_val = np.average(np.diff(self.depths_layer))

        try:
            self.model_array = _search_vel_file(velmod)
        except (IOError,ValueError):
            raise IOError
        else:
            if layer_mod:
                self.vp, self.vs, self.rho = \
                    _layer2grid(self.depths_layer, self.model_array)
            else:
                self.vp, self.vs, self.rho = \
                    _discretize(self.depths_layer, self.model_array)


    @classmethod
    def read_layer_model(cls, dep_range, h, vp, vs, rho=None, elevation=0):
        mod = cls(dep_range, velmod=None, layer_mod=True, elevation=elevation)
        return mod


    def plot_model(self, show=True):
        plt.style.use('bmh')
        if self.isrho:
            self.model_fig = plt.figure(figsize=(6,6))
            fignum = 2
        else:
            self.model_fig = plt.figure(figsize=(4,6))
            fignum = 1
        self.model_ax = self.model_fig.add_subplot(1,fignum,1)
        self.model_ax.step(self.vp, self.depths_layer, where='pre', label='Vp')
        self.model_ax.step(self.vs, self.depths_layer, where='pre', label='Vs')
        self.model_ax.legend()
        self.model_ax.set_xlabel('Velocity (km/s)')
        self.model_ax.set_ylabel('Depth (km)')
        self.model_ax.set_ylim([self.depths_layer[0], self.depths_layer[-1]])
        self.model_ax.invert_yaxis()
        if self.isrho:
            self.rho_ax = self.model_fig.add_subplot(1,fignum,2)
            self.rho_ax.step(self.rho, self.depths_layer, where='pre', color='C2', label='Density')
            self.rho_ax.legend()
            self.rho_ax.set_xlabel('Density (km/s)')
            self.rho_ax.set_ylim([self.depths_layer[0], self.depths_layer[-1]])
            self.rho_ax.invert_yaxis()
        if show:
            plt.show()




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
        """
        calculate piercing point, P for Sp and S for Ps
        Parameters
        ----------
        rayp
        phase
        sphere

        Returns
        -------

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

if __name__ == "__main__":
    depth = np.array([0, 20.1, 35.1, 100])
    model = DepModel(depth)
    rayp = np.arange(0.04, 0.09, 0.01)
    print(model.vp)
