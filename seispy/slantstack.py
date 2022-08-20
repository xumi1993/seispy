import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from seispy.geo import srad2skm, skm2sdeg

class SlantStack():
    def __init__(self, seis, timeaxis, dis) -> None:
        """Slant stack in tau domain. Refer to Tauzin et al., 2008 JGR in detail.

        :param seis: 2D array for RFs with shape of (ev_num, npts)
        :type seis: numpy.ndarray
        :param timeaxis: 1D array for time axis
        :type timeaxis: numpy.ndarray
        :param dis: 1D array for all event distance in km
        :type dis: numpy.ndarray
        """
        self.datar = seis
        self.time_axis = timeaxis
        self.dis = dis
        self.ref_dis = 65
        self.rayp_range = np.arange(-0.35, 0.35, 0.01)
        self.tau_range = np.arange(0, 100, 0.1)
        self.syn_tau = np.array([])
        self.syn_drayp = np.array([])

    def stack(self, ref_dis=None, rayp_range=None, tau_range=None):
        if ref_dis is not None and isinstance(ref_dis, (int, float)):
            self.ref_dis = ref_dis
        elif ref_dis is None:
            pass
        else:
            raise TypeError('{} should be in int or float type.'.format(ref_dis))
        if rayp_range is not None and isinstance(rayp_range, np.ndarray):
            self.rayp_range = rayp_range
        elif rayp_range is None:
            pass
        else:
            raise TypeError('{} should be in numpy.ndarray type.'.format(rayp_range))
        if tau_range is not None and isinstance(tau_range, np.ndarray):
            self.tau_range = tau_range
        elif tau_range is None:
            pass
        else:
            raise TypeError('{} should be in numpy.ndarray type.'.format(tau_range))
        ev_num = self.datar.shape[0]
        taus, rayps = np.meshgrid(self.tau_range, self.rayp_range)
        self.stack_amp = np.zeros([self.rayp_range.shape[0], self.tau_range.shape[0]])
        for i in range(ev_num):
            tps = taus - rayps * (self.dis[i] - self.ref_dis)
            self.stack_amp += interp1d(self.time_axis, self.datar[i, :], fill_value='extrapolate')(tps)
        self.stack_amp /= ev_num

    def syn_tps(self, phase_list, velmodel='iasp91', focal_dep=10):
        model = TauPyModel(model=velmodel)
        phase_list.insert(0, 'P')
        arrs = model.get_travel_times(focal_dep, self.ref_dis, phase_list=phase_list)
        p_arr = arrs[0].time
        p_rayp = skm2sdeg(srad2skm(arrs[0].ray_param))
        self.syn_tau = [arr.time - p_arr for arr in arrs[1:]]
        self.syn_drayp = [p_rayp - skm2sdeg(srad2skm(arr.ray_param))for arr in arrs[1:]]

    def plot(self, cmap='jet', xlim=None, vmin=None, vmax=None, figpath=None, colorbar=True):
        """ Imaging for slant stacking

        Parameters
        ----------
        cmap : str, optional
            The colormap to use, by default 'jet'
        xlim : _type_, optional
             Set the x limits of the current axes., by default None
        vmin : _type_, optional
            Normalize the minium value data to the ``vmin``, by default None
        vmax : _type_, optional
            Normalize the maximum value data to the ``vmax``, by default None
        figpath : _type_, optional
            Output path to the image, by default None
        colorbar : bool, optional
            Wether plot the colorbar, by default True
        """
        plt.style.use("bmh")
        self.fig = plt.figure(figsize=(8,5))
        self.ax = self.fig.add_subplot()
        if vmin is None and vmax is None:
            vmax = np.max(np.abs(self.stack_amp))*0.1
            vmin = -vmax
        elif vmin is None and isinstance(vmax, (int, float)):
            vmin = np.min(self.stack_amp)
        elif vmax is None and isinstance(vmin, (int, float)):
            vmax = np.max(self.stack_amp)
        elif isinstance(vmax, (int, float)) and isinstance(vmin, (int, float)):
            pass
        else:
            raise TypeError('vmin and vmax must be in int or in float type')
        im = self.ax.pcolor(self.tau_range, self.rayp_range, self.stack_amp,
                            vmax=vmax, vmin=vmin, cmap=cmap)
        if colorbar:
            cax = self.fig.colorbar(im, ax=self.ax)
            cax.set_label('Stack Amplitude')
        self.ax.grid()
        self.ax.scatter(self.syn_tau, self.syn_drayp, color='k', marker='x')
        if xlim is not None and isinstance(xlim, (list, np.ndarray)):
            self.ax.set_xlim(xlim)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Slowness (s/$^\circ$)')
        if figpath is not None and isinstance(figpath, str):
            self.fig.savefig(figpath, format='png', dpi=400, bbox_inches='tight')
        else:
            plt.show()
