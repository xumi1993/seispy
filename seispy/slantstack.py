import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from seispy.geo import srad2skm, skm2sdeg
from seispy.utils import scalar_instance, array_instance

class SlantStack():
    def __init__(self, seis, timeaxis, dis) -> None:
        """Slant stack in tau domain. Refer to Tauzin et al., 2008 JGR in detail.

        :param seis: 2D array for RFs with shape of (ev_num, npts)
        :type seis: numpy.ndarray
        :param timeaxis: 1D array for time axis
        :type timeaxis: numpy.ndarray
        :param dis: 1D array for all event distance in deg
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
        """Slant stack for receiver function

        :param ref_dis: reference distance, by default None
        :type ref_dis: int or float, optional
        :param rayp_range: range of ray parameter, by default None
        :type rayp_range: numpy.ndarray, optional
        :param tau_range: range of tau, by default None
        :type tau_range: numpy.ndarray, optional
        """
        if ref_dis is not None and scalar_instance(ref_dis):
            self.ref_dis = ref_dis
        elif ref_dis is None:
            pass
        else:
            raise TypeError('{} should be in int or float type.'.format(ref_dis))
        if rayp_range is not None and array_instance(rayp_range):
            self.rayp_range = rayp_range
        elif rayp_range is None:
            pass
        else:
            raise TypeError('{} should be in numpy.ndarray type.'.format(rayp_range))
        if tau_range is not None and array_instance(tau_range):
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
        """Calculate the theoretical tau and reference rayp for the given phase list

        :param phase_list: phase list for theoretical arrival time
        :type phase_list: list
        :param velmodel: velocity model, by default 'iasp91'
        :type velmodel: str, optional
        :param focal_dep: focal depth, by default 10
        :type focal_dep: int or float, optional
        """
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
        elif vmin is None and scalar_instance(vmax):
            vmin = np.min(self.stack_amp)
        elif vmax is None and scalar_instance(vmin):
            vmax = np.max(self.stack_amp)
        elif scalar_instance(vmax) and scalar_instance(vmin):
            pass
        else:
            raise TypeError('vmin and vmax must be in int or in float type')
        im = self.ax.pcolor(self.tau_range, self.rayp_range, self.stack_amp,
                            vmax=vmax, vmin=vmin, cmap=cmap)
        if colorbar:
            cax = self.fig.colorbar(im, ax=self.ax)
            cax.set_label('Stack Amplitude')
        self.ax.grid(visible=True)
        self.ax.scatter(self.syn_tau, self.syn_drayp, color='k', marker='x')
        if xlim is not None and array_instance(xlim):
            self.ax.set_xlim(xlim)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Slowness (s/$^\circ$)')
        if figpath is not None and isinstance(figpath, str):
            self.fig.savefig(figpath, format='png', dpi=400, bbox_inches='tight')
        else:
            plt.show()
