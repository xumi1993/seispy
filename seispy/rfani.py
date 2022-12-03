import numpy as np
from obspy.io.sac import SACTrace
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
from os.path import join, abspath, dirname
from seispy.geo import cosd, sind, extrema
from scipy.interpolate import griddata
from seispy.utils import load_cyan_map


def joint_stack(energy_r, energy_cc, energy_tc, weight=[0.4, 0.3, 0.3]):
    energy_r = energy_r / np.max(energy_r)
    energy_cc = energy_cc / np.max(energy_cc)
    energy_tc = energy_tc / np.max(energy_tc)
    return np.exp(np.log(energy_r) * weight[0] + np.log(energy_cc) * weight[1] - np.log(energy_tc) * weight[2])


def average_delay(fd, td):
    uniq_fd = np.unique(fd)
    if uniq_fd.size > 2:
        raise ValueError('FVD do not converge, {} gotten'.format(uniq_fd))
    uniq_td = []
    for i, fv in enumerate(uniq_fd):
        idx = np.where(fd==fv)[0]
        uniq_td.append(np.mean(td[idx]))
    return uniq_fd, np.array(uniq_td)


class RFAni():
    def __init__(self, sacdatar, tb, te, tlen=3, val=10, rayp=0.06, model='iasp91'):
        """ Estimate crustal anisotropy with a joint method. 
            See Liu and Niu (2012, doi: 10.1111/j.1365-246X.2011.05249.x) in detail.

        Parameters
        ----------
        sacdatar : seispy.rfcorrect.RFStation
            RFStation include RF data
        tb : float
            Starting time for searching Ps peak.
        te : float
            end time for searching Ps peak.
        tlen : int, optional
            Half-length of time window for trimming waveform around Ps, by default 3
        val : int, optional
            Interval of Back-azimuth, by default 10
        rayp : float, optional
            Reference ray-parameter for moveout correction, by default 0.06
        model : str, optional
            1D velocity for moveout correction, by default 'iasp91'
        """
        self.tb = tb
        self.te = te
        self.tlen = tlen
        self.sacdatar = sacdatar
        self.sacdatar.resample(0.1)
        self.nbs = int((self.tb + self.sacdatar.shift) / self.sacdatar.sampling)
        self.nes = int((self.te + self.sacdatar.shift) / self.sacdatar.sampling)
        self.sacdatar.moveoutcorrect(ref_rayp=rayp, velmod=model, replace=True)
        self.baz_stack(val=val)
        self.search_peak_amp()
        self.init_ani_para()
        self.fvd, self.deltat = np.meshgrid(self.fvd_1d, self.deltat_1d)

    def baz_stack(self, val=10):
        self.stack_range = np.arange(0, 360, val)
        self.rft_baz = np.zeros([self.stack_range.shape[0], self.sacdatar.rflength])
        self.rfr_baz = np.zeros_like(self.rft_baz)
        self.count_baz = np.zeros(self.stack_range.shape[0])
        search_range = np.append(self.stack_range, self.stack_range[-1]+val)
        for i in range(self.stack_range.size):
            idx = np.where((self.sacdatar.bazi > search_range[i]) & (self.sacdatar.bazi < search_range[i+1]))[0]
            self.count_baz[i] = idx.size
            if idx.size != 0:
                self.rft_baz[i] = np.mean(self.sacdatar.datat[idx], axis=0)
                self.rfr_baz[i] = np.mean(self.sacdatar.datar[idx], axis=0)

    def search_peak_amp(self):
        mean_rf = np.mean(self.rfr_baz, axis=0)
        nmax = extrema(mean_rf[self.nbs:self.nes])+self.nbs
        if nmax.size > 1:
            nps = nmax[np.nanargmax(mean_rf[nmax])]
        else:
            nps = nmax[0]
        self.nb = int(nps - self.tlen / self.sacdatar.sampling)
        self.ne = int(nps + self.tlen / self.sacdatar.sampling)

    def init_ani_para(self):
        self.deltat_1d = np.arange(0, 1.55, 0.05)
        self.fvd_1d = np.arange(0, 365, 5)

    def cut_energy_waveform(self, idx, nb, ne):
        engr = np.zeros([nb.shape[0], nb.shape[1], self.ne-self.nb])
        for i in range(nb.shape[0]):
            for j in range(nb.shape[1]):
                engr[i, j, :] = self.rfr_baz[idx, nb[i, j]:ne[i, j]]
        return engr

    def radial_energy_max(self):
        energy = np.zeros([self.fvd.shape[0], self.fvd.shape[1], self.ne-self.nb])
        # tmp_data = np.zeros(self.ne-self.nb)
        for i, baz in enumerate(self.stack_range):
            t_corr = (self.deltat / 2) * cosd(2 * (self.fvd - baz))
            nt_corr = (t_corr / self.sacdatar.sampling).astype(int)
            new_nb = self.nb - nt_corr
            new_ne = self.ne - nt_corr
            energy += self.cut_energy_waveform(i, new_nb, new_ne)
        energy = np.max(energy ** 2, axis=2)
        energy /= np.max(np.sum(self.rfr_baz[:, self.nb:self.ne], axis=0)**2)
        return energy

    def xyz2grd(self, energy):
        self.fvd, self.deltat = np.meshgrid(self.fvd_1d, self.deltat_1d)
        return griddata(self.ani_points, energy, (self.fvd, self.deltat))

    def rotate_to_fast_slow(self):
        self.ani_points = np.empty([0, 2])
        for f in self.fvd_1d:
            for d in self.deltat_1d:
                self.ani_points = np.vstack((self.ani_points, np.array([f, d])))
        energy_cc = np.zeros(self.ani_points.shape[0])
        energy_tc = np.zeros(self.ani_points.shape[0])
        raw_energy_r = np.sum(np.sum(self.rfr_baz[:, self.nb:self.ne], axis=0) ** 2 - np.sum(self.rfr_baz[:, self.nb:self.ne] ** 2, axis=0))
        raw_energy_t = np.sum(np.sum(self.rft_baz[:, self.nb:self.ne] ** 2, axis=0))
        for i, point in enumerate(self.ani_points):
            nt_corr = (point[1]/2 / self.sacdatar.sampling).astype(int)
            nt_fast = np.arange(self.nb, self.ne) + nt_corr
            nt_slow = np.arange(self.nb, self.ne) - nt_corr
            fcr = np.zeros(self.ne-self.nb)
            fcr_sq = np.zeros(self.ne-self.nb)
            fct = 0
            for j, baz in enumerate(self.stack_range):
                data_fast = self.rfr_baz[j, nt_slow] * cosd(point[0] - baz) + self.rft_baz[j, nt_slow] * sind(point[0] - baz)
                data_slow = -self.rfr_baz[j, nt_fast] * sind(point[0] - baz) + self.rft_baz[j, nt_fast] * cosd(point[0] - baz)
                back_rotate_data_r = data_fast * cosd(point[0] - baz) - data_slow * sind(point[0] - baz)
                back_rotate_data_t = data_fast * sind(point[0] - baz) + data_slow * cosd(point[0] - baz)
                fcr += back_rotate_data_r
                fcr_sq += back_rotate_data_r ** 2
                # print(np.sum(back_rotate_data_t ** 2), np.sum(self.rft_baz[j, self.nb:self.ne] ** 2))
                fct += np.sum(back_rotate_data_t ** 2)
            energy_cc[i] = np.sum(fcr ** 2 - fcr_sq) / raw_energy_r
            energy_tc[i] = fct/ raw_energy_t
        # energy_cc /= np.max(np.abs(energy_cc))
        energy_tc /= np.max(np.abs(energy_tc))
        return self.xyz2grd(energy_cc), self.xyz2grd(energy_tc)

    def plot_stack_baz(self, enf=60, outpath='./'):
        ml = MultipleLocator(5)
        bound = np.zeros_like(self.sacdatar.time_axis)
        plt.style.use("bmh")
        plt.rc('grid', color='white', linestyle='-', linewidth=0.7)
        plt.rcParams["axes.grid.axis"] = "x"
        fig = plt.figure(figsize=(15, 8))
        axr = plt.subplot(1, 2, 1)
        for i, baz in enumerate(self.stack_range):
            if np.mean(self.rfr_baz[i]) == 0:
                continue
            amp = self.rfr_baz[i] * enf + baz
            # axr.plot(time_axis, amp, linewidth=0.2, color='black')
            axr.fill_between(self.sacdatar.time_axis, amp, bound + baz, where=amp > self.stack_range[i], facecolor='red', alpha=0.7)
            axr.fill_between(self.sacdatar.time_axis, amp, bound + baz, where=amp < self.stack_range[i], facecolor='#1193F4', alpha=0.7)
        # axr.plot([0, 0], [0, 360], linewidth=1, color='black')
        # axr.set_xlim([-1, 15])
        # axr.xaxis.set_major_locator(MultipleLocator(2))
        # axr.set_ylim(0, 360)
        # axr.set_yticks(np.arange(0, 360+30, 30))
        # axr.yaxis.set_minor_locator(ml)
        # axr.set_ylabel('Back-azimuth ($^\circ$)')
        # axr.set_xlabel('Time after P(s)')
        axr.set_title('{} component ({})'.format(self.sacdatar.comp, self.sacdatar.staname), fontsize=16)

        axt = plt.subplot(1, 2, 2)
        # axt.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        for i, baz in enumerate(self.stack_range):
            if np.mean(self.rft_baz[i]) == 0:
                continue
            amp = self.rft_baz[i] * enf + baz
            # axt.plot(time_axis, amp, linewidth=0.2, color='black')
            axt.fill_between(self.sacdatar.time_axis, amp, bound + baz, where=amp > baz, facecolor='red', alpha=0.7)
            axt.fill_between(self.sacdatar.time_axis, amp, bound + baz, where=amp < baz, facecolor='#1193F4', alpha=0.7)
        for ax in [axr, axt]:
            ax.plot([0, 0], [-30, 400], linewidth=1, color='black')
            ax.set_xlim(-1, 15)
            ax.xaxis.set_major_locator(MultipleLocator(2))
            ax.set_ylim(-10, 370)
            ax.set_yticks(np.arange(0, 360+30, 30))
            ax.yaxis.set_minor_locator(ml)
            ax.set_xlabel('Time after P (s)',  fontsize=16)
        # axt.plot([0, 0], [0, 360], linewidth=1, color='black')
        # axt.set_xlim([-1, 15])
        # # axt.set_xticks(np.arange(0, 5))
        # axt.xaxis.set_major_locator(MultipleLocator(2))
        # axt.set_ylim(0, 360)
        # axt.set_yticks(np.arange(0, 360+30, 30))
        # axt.yaxis.set_minor_locator(ml)
        axr.set_ylabel('Back-azimuth ($^\circ$)',  fontsize=16)
        # axt.set_xlabel('Time after P(s)')
        axt.set_title('T component ({})'.format(self.sacdatar.staname), fontsize=16)
        fig.savefig(join(outpath, '{}_baz_stack.png'.format(self.sacdatar.staname)), dpi=400, bbox_inches='tight')

    def plot_correct(self, fvd=0, dt=0.44, enf=80, outpath=None):
        nt_corr = int((dt/2 / self.sacdatar.sampling))
        # nt_fast = np.arange(self.nb, self.ne) + nt_corr
        # nt_slow = np.arange(self.nb, self.ne) - nt_corr
        time_axis = np.arange(self.nb, self.ne) * self.sacdatar.sampling - self.sacdatar.shift
        bound = np.zeros_like(time_axis)
        ml = MultipleLocator(5)
        plt.figure(figsize=(8, 6))
        axr = plt.subplot(1, 2, 1)
        axt = plt.subplot(1, 2, 2)
        for j, baz in enumerate(self.sacdatar.baz):
            rot_fast = self.rfr_baz[j] * cosd(fvd - baz) + self.rft_baz[j] * sind(fvd - baz)
            rot_slow = -self.rfr_baz[j] * sind(fvd - baz) + self.rft_baz[j] * cosd(fvd - baz)
            data_fast = rot_fast[self.nb-nt_corr:self.ne-nt_corr]
            data_slow = rot_slow[self.nb+nt_corr:self.ne+nt_corr]
            back_rotate_data_r = data_fast * cosd(fvd - baz) - data_slow * sind(fvd - baz)
            back_rotate_data_t = data_fast * sind(fvd - baz) + data_slow * cosd(fvd - baz)
            amp = back_rotate_data_r * enf + baz
            axr.plot(time_axis, amp, linewidth=0.2, color='black')
            axr.fill_between(time_axis, amp, bound + baz, where=amp > baz, facecolor='red', alpha=0.5)
            axr.fill_between(time_axis, amp, bound + baz, where=amp < baz, facecolor='blue', alpha=0.5)
            amp = back_rotate_data_t * enf + baz
            axt.plot(time_axis, amp, linewidth=0.2, color='black')
            axt.fill_between(time_axis, amp, bound + baz, where=amp > baz, facecolor='red', alpha=0.5)
            axt.fill_between(time_axis, amp, bound + baz, where=amp < baz, facecolor='blue', alpha=0.5)
        for ax in [axr, axt]:
            ax.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
            ax.set_xlim(-1, 10)
            ax.set_xticks(np.arange(0, 11, 1))
            ax.set_ylim(-20, 370)
            ax.set_yticks(np.arange(0, 360+30, 30))
            ax.yaxis.set_minor_locator(ml)
            ax.set_xlabel('Time after P(s)',  fontsize=16)
        axr.set_ylabel('Back-azimuth ($^\circ$)', fontsize=16)
        axr.set_title('R component',  fontsize=16)
        axt.set_title('T component',  fontsize=16)
        if outpath is not None and isinstance(outpath, str):
            plt.savefig(join(outpath, 'rf_corrected.png'), dpi=400, bbox_inches='tight')

    def search_peak_list(self, energy, opt='max'):
        if opt == 'max':
            ind = np.argwhere(energy == np.max(energy))
        elif opt == 'min':
            ind = np.argwhere(energy == np.min(energy))
        else:
            raise ValueError('\'opt\' must be max or min')
        return self.ani_points[ind][:, 0], self.ani_points[ind][:, 1]

    def search_peak(self, energy, opt='max'):
        if opt == 'max':
            ind = np.argwhere(energy == np.max(energy))
        elif opt == 'min':
            ind = np.argwhere(energy == np.min(energy))
        else:
            raise ValueError('\'opt\' must be max or min')
        best_fvd = []
        best_dt = []
        for i, j in ind:
            best_fvd.append(self.fvd[i, j])
            best_dt.append(self.deltat[i, j])
        uniq_fd, uniq_td = average_delay(np.array(best_fvd), np.array(best_dt))
        return uniq_fd, uniq_td

    def joint_ani(self, weight=[0.5, 0.3, 0.2]):
        self.energy_r = self.radial_energy_max()
        self.energy_cc, self.energy_tc = self.rotate_to_fast_slow()
        self.energy_joint = joint_stack(self.energy_r, self.energy_cc, self.energy_tc, weight)
        self.bf, self.bt = self.search_peak(self.energy_joint, opt='max')
        return self.bf, self.bt

    def plot_polar(self, cmap=load_cyan_map(), show=False, outpath='./'):
        """Polar map of crustal anisotropy. See Liu and Niu (2012, doi: 10.1111/j.1365-246X.2011.05249.x) in detail.

        :param cmap: Colormap of matplotlib, defaults to 'rainbow'
        :type cmap: optional
        :param show: If show the polar map in the Matplotlib window, defaults to True
        :type show: bool, optional
        :param outpath: Output path to saving the figure. If show the figure in the Matplotlib window, this option will be invalid, defaults to current directory.
        :type outpath: str, optional
        """
        fig, axes = plt.subplots(2, 2, figsize=(8, 7), subplot_kw={'projection': 'polar'}, constrained_layout=True)
        axs = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]
        energy_all = [self.energy_r, self.energy_cc, self.energy_tc, self.energy_joint]
        energy_title = ['R cosine energy', 'R cross-correlation', 'T energy', 'Joint']
        for ax, energy, title in zip(axs, energy_all, energy_title):
            ax.set_theta_direction(-1)
            ax.set_theta_zero_location("N")
            if title == 'T energy':
                eng = ax.pcolor(np.radians(self.fvd), self.deltat, energy, cmap=cmap.reversed(), shading='auto')
            else:
                eng = ax.pcolor(np.radians(self.fvd), self.deltat, energy, cmap=cmap, shading='auto')
            ax.grid(True, color='lightgray', linewidth=0.5)
            ax.scatter(np.radians(self.bf), self.bt, color='white', marker='X', s=48)
            ax.set_xticks(np.radians(np.arange(0, 360, 30)))
            ax.set_yticks(np.arange(0, 1.5, 0.5))
            ax.set_title(title)
            fig.colorbar(eng, ax=ax)
        fig.suptitle('{}\nFVD = ${:.0f}^\circ$, $\delta t$ = ${:.1f}$ s'.format(self.sacdatar.staname, self.bf[0], self.bt[0]), fontsize=16)
        if show:
            plt.show()
        else:
            fig.savefig(join(outpath, self.sacdatar.staname+'_joint_ani.png'), dpi=400, bbox_inches='tight')


if __name__ == "__main__":
    pass
