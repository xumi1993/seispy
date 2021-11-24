import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

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
        self.tau_range = np.arange(0, 100)

    def stack(self, ref_dis=None, rayp_range=None, tau_range=None):
        if ref_dis is not None and isinstance(ref_dis, (int, float)):
            self.ref_dis = ref_dis
        else:
            raise TypeError('{} should be in int or float type.'.format(ref_dis))
        if rayp_range is not None and isinstance(rayp_range, np.ndarray):
            self.rayp_range = rayp_range
        else:
            raise TypeError('{} should be in numpy.ndarray type.'.format(rayp_range))
        if tau_range is not None and isinstance(tau_range, np.ndarray):
            self.tau_range = tau_range
        else:
            raise TypeError('{} should be in numpy.ndarray type.'.format(tau_range))
        ev_num = self.datar.shape[0]
        tmp = np.zeros([ev_num, tau_range.shape[0]])
        self.stack_amp = np.zeros([rayp_range.shape[0], tau_range.shape[0]])
        for j in range(rayp_range.shape[0]):
            for i in range(ev_num):
                self.datar[i, :] = self.datar[i, :] / np.max(np.abs(self.datar[i, :]))
                tps = tau_range - rayp_range[j] * (self.dis[i] - ref_dis)
                tmp[i, :] = interp1d(self.time_axis, self.datar[i, :], fill_value='extrapolate')(tps)
            self.stack_amp[j, :] = np.mean(tmp, axis=0)
        
    def plot(self, cmap='jet', xlim=None, vmin=None, vmax=None, figpath=None):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()
        if vmin is None and vmax is None:
            vmax = np.max(np.abs(self.stack_amp))
            vmin = -vmax
        elif vmin is None and isinstance(vmax, (int, float)):
            vmin = np.min(self.stack_amp)
        elif vmax is None and isinstance(vmin, (int, float)):
            vmax = np.max(self.stack_amp)
        elif isinstance(vmax, (int, float)) and isinstance(vmin, (int, float)):
            pass
        else:
            raise TypeError('vmin and vmax must be in int or in float type')
        cs = self.ax.pcolor(self.tau_range, self.rayp_range, self.stack_amp,
                            vmax=vmax, vmin=vmin, cmap=cmap)
        self.ax.colorbar(cs)
        if xlim is not None and isinstance(xlim, (list, np.ndarray)):
            self.ax.set_xlim(xlim)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Slowness (s/km)')
        if figpath is not None and isinstance(figpath, str):
            self.fig.savefig(figpath, format='png', dpi=400, bbox_inches='tight')
        else:
            plt.show()
