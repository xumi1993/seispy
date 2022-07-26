import numpy as np
from scipy.sparse.linalg import lsqr
from seispy.geo import cosd, sind
import matplotlib.pyplot as plt
from obspy.io.sac import SACTrace
from os.path import join



class Harmonics():
    def __init__(self, rfsta, tmin=-5, tmax=10) -> None:
        """ Harmonic decomposition for extracting anisotropic and isotropic features from the radial and transverse RFs

        :param rfsta: data class of RFStation
        :type rfsta: :class:`seispy.rfcorrect.RFStation`
        :param tmin: Start time relative to P
        :type tmin: float
        :param tmax: End time relative to P
        :type tmax: float
        """
        self.rfsta = rfsta
        self.tmin = tmin
        self.tmax = tmax
        self.cut_trace()

    def cut_trace(self):
        """Trim RFs from tmin to tmax
        """
        nb = int((self.tmin + self.rfsta.shift) / self.rfsta.sampling)
        ne = int((self.tmax + self.rfsta.shift) / self.rfsta.sampling)
        self.datar = self.rfsta.datar[:, nb:ne+1]
        self.datat = self.rfsta.datat[:, nb:ne+1]
        self.nsamp = ne - nb + 1
        self.time_axis = np.linspace(self.tmin, self.tmax, self.nsamp)

    def harmo_trans(self):
        self.traces = np.vstack((self.datar, self.datat))
        harmonic = np.zeros((self.rfsta.ev_num*2, 5))
        unmodel = np.zeros((self.rfsta.ev_num*2, 5))
        self.harmonic_trans = np.zeros((5, self.nsamp))
        self.unmodel_trans = np.zeros((5, self.nsamp))
        for i, baz in enumerate(self.rfsta.bazi):
            harmonic[i] = np.array([1, cosd(baz), sind(baz), cosd(baz*2), sind(baz*2)])
            harmonic[i+self.rfsta.ev_num] = np.array([0, cosd(baz+90), sind(baz+90), cosd(baz*2+90), sind(baz*2+90)])
            unmodel[i] = np.array([1, cosd(baz), sind(baz), cosd(baz*2), sind(baz*2)])
            unmodel[i+self.rfsta.ev_num] = np.array([0, cosd(baz-90), sind(baz-90), cosd(baz*2-90), sind(baz*2-90)])
        for i in range(self.nsamp):
            b = self.traces[:, i]
            self.harmonic_trans[:, i] = lsqr(harmonic, b, damp=0)[0]
            self.unmodel_trans[:, i] = lsqr(unmodel, b, damp=0)[0]
    
    def write_constant(self, out_sac_path='./'):
        """ Write constant component to SAC file.

        :param out_sac_path: Output path, defaults to './'
        :type out_sac_path: str, optional
        """
        sac = SACTrace(data=self.harmonic_trans[0, :])
        sac.b = self.tmin
        sac.delta = self.rfsta.sampling
        sac.user0 = 0.0
        sac.user1 = self.rfsta.f0[0]
        sac.stla = self.rfsta.stla
        sac.stlo = self.rfsta.stlo
        sac.stel = self.rfsta.stel
        sac.write(join(out_sac_path, '{}_constant_R.sac'.format(self.rfsta.staname)))

    def plot(self, outpath='./', enf=2.):
        """Plot harmonic and unmodeled components

        :param outpath: Output path, defaults to './'
        :type outpath: str, optional
        :param enf: Amplification factor, defaults to 2.0
        :type enf: float, optional
        """
        plt.style.use("bmh")
        plt.rc('grid', color='white', linestyle='-', linewidth=0.7)
        plt.rcParams["axes.grid.axis"] = "x"
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        mtx = [self.harmonic_trans, self.unmodel_trans]
        bound = np.zeros(self.nsamp)
        titles = ['Dipping/Anisotropic', 'Unmodeled']
        for iax, ax in enumerate(axes):
            for i in range(5):
                data = mtx[iax][5-i-1] * enf + (i + 1)
                ax.fill_between(self.time_axis, data, bound + i+1, where=data > i+1, facecolor='red',
                                alpha=0.7)
                ax.fill_between(self.time_axis, data, bound + i+1, where=data < i+1, facecolor='#1193F4',
                                alpha=0.7)
            ax.plot([0, 0], [0, 6], color='k', lw=1.2)
            ax.set_xlim([-1, self.tmax])
            ax.set_xlabel('Time after P (s)', fontsize=13)
            ax.set_ylim([0, 6])
            ax.set_title(titles[iax])
        axes[0].set_yticks([1, 2, 3, 4, 5])
        axes[0].set_yticklabels(['sin2${\\theta}$', 'cos2${\\theta}$', 'sin${\\theta}$', 'cos${\\theta}$', 'Constant'], fontsize=13)
        fig.savefig(join(outpath, '{}_harmonic_trans.png'.format(self.rfsta.staname)), dpi=300, bbox_inches='tight')
