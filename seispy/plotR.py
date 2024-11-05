import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
from os.path import join


class PlotR:
    def __init__(self, stadata, enf=12, xlim=[2, 80], key='bazi'):
        """Plot receiver function in R component

        :param stadata: receiver function data
        :type stadata: sespy.rfcorrect.RFStation
        :param enf: enlarge factor, defaults to 12
        :type enf: int, optional
        :param xlim: xlim, defaults to [2, 80]
        :type xlim: list, optional
        :param key: sort key, defaults to 'bazi'
        :type key: str, optional
        """
        self.stadata = stadata
        self.stadata.sort(key)
        self.enf = enf
        self.xlim = xlim
        self.key = key
        self.init_figure()

    def init_figure(self):
        """Initialize figure
        """
        self.fig = plt.figure(figsize=(8, 10))
        self.fig.tight_layout()
        gs = GridSpec(15, 3)
        gs.update(wspace=0.25)
        self.axr = plt.subplot(gs[1:, 0:-1])
        self.axr.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        self.axb = plt.subplot(gs[1:, -1])
        self.axb.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        self.axs = plt.subplot(gs[0, 0:-1])
        self.axs.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')

    def plot_waves(self):
        """Plot PRFs with R component
        """
        bound = np.zeros(self.stadata.rflength)
        for i in range(self.stadata.ev_num):
            datar = self.stadata.data_prime[i] * self.enf + (i + 1)
            self.axr.fill_between(self.stadata.time_axis, datar, bound + i+1, where=datar > i+1, facecolor='red',
                            alpha=0.7)
            self.axr.fill_between(self.stadata.time_axis, datar, bound + i+1, where=datar < i+1, facecolor='blue',
                            alpha=0.7)
        # rfcorr, _ = stadata.moveoutcorrect( 0.1, np.arange(300), sphere=False)
        # rfsum = np.mean(rfcorr, axis=0)
        rfsum = np.mean(self.stadata.data_prime, axis=0)
        self.axs.fill_between(self.stadata.time_axis, rfsum, 0, where=rfsum > 0, facecolor='red', alpha=0.7)
        self.axs.fill_between(self.stadata.time_axis, rfsum, 0, where=rfsum < 0, facecolor='blue', alpha=0.7)
        # axs.plot(stadata.time_axis, rfsum, color='gray', lw=0.5)
        self.axb.scatter(self.stadata.bazi, np.arange(self.stadata.ev_num) + 1, s=7)
        # axp = axb.twiny()
        # axp.scatter(stadata.rayp, np.arange(stadata.ev_num) + 1, s=7)
        # return axp

    def set_fig(self):
        """Set figure
        """
        y_range = np.arange(self.stadata.ev_num) + 1
        x_range = np.arange(0, self.xlim[1]+2, 5)
        space = 2

        # set axr
        self.axr.set_xlim(self.xlim)
        self.axr.set_xticks(x_range)
        self.axr.set_xticklabels(x_range, fontsize=8)
        self.axr.set_ylim(0, self.stadata.ev_num + space)
        self.axr.set_yticks(y_range)
        self.axr.set_yticklabels(self.stadata.event, fontsize=5)
        self.axr.set_xlabel('Time after P (s)', fontsize=13)
        self.axr.set_ylabel('Event', fontsize=13)
        self.axr.add_line(Line2D([0, 0], self.axr.get_ylim(), color='black'))
        # axr.set_title('R components ({})'.format(stadata.staname), fontsize=16)

        # set axb
        self.axb.set_xlim(0, 360)
        self.axb.set_xticks(np.linspace(0, 360, 7))
        self.axb.set_xticklabels(np.linspace(0, 360, 7, dtype='i'), fontsize=8)
        self.axb.set_ylim(0, self.stadata.ev_num + space)
        self.axb.set_yticks(y_range)
        self.axb.set_yticklabels(y_range, fontsize=5)
        self.axb.set_xlabel(r'Back-azimuth ($^\circ$)', fontsize=13)
        # axp.set_xlabel('Ray-parameter (s/km)', fontsize=13)

        self.axs.set_title('{} components ({})'.format(self.stadata.comp.upper(), self.stadata.staname), fontsize=16)
        self.axs.set_xlim(self.xlim)
        self.axs.set_xticks(x_range)
        self.axs.set_xticklabels([])
        self.axs.set_yticks([np.sum(self.axs.get_ylim())/3])
        self.axs.tick_params(axis='y', left=False)
        self.axs.set_yticklabels(['Sum'], fontsize=8)
        self.axs.add_line(Line2D([0, 0], self.axs.get_ylim(), color='black'))

    def plot(self, out_path=None, outformat='g', show=False):
        """ Plot receiver function in R component
        :param out_path: output path
        :type out_path: str
        :param outformat: output format
        :type outformat: str
        :param show: show figure
        :type show: bool
        """
        self.plot_waves()
        self.set_fig()
        if out_path is not None:
            if outformat == 'g':
                self.fig.savefig(join(out_path, self.stadata.staname+'_{}_{}order_{:.1f}.png'.format(
                                self.stadata.comp, self.key, self.stadata.f0[0])),
                                dpi=400, bbox_inches='tight')
            elif outformat == 'f':
                self.fig.savefig(join(out_path, self.stadata.staname+'_{}_{}order_{:.1f}.pdf'.format(
                                self.stadata.comp, self.key, self.stadata.f0[0])),
                                format='pdf', bbox_inches='tight')
        if show:
            plt.show()


def plotr(rfsta, out_path='./', xlim=[-2, 80], key='bazi', enf=6, outformat='g', show=False):
    """ Plot receiver function in R component
    :param rfsta: Path to PRFs
    :type rfsta: seispy.rfcorrect.RFStation
    :param enf: The enlarge factor, defaults to 6
    :type enf: int, optional
    :param out_path: The output path, defaults to current directory
    :type out_path: str, optional
    :param key: The key to sort PRFs, available for ``event``, ``evla``, ``evlo``, ``evdp``,
        ``dis``, ``bazi``, ``rayp``, ``mag``, ``f0``, defaults to ``bazi``
    :type key: str, optional
    :param outformat: File format of the image file, g as \'png\', f as \'pdf\', defaults to 'g'
    :type outformat: str, optional
    :param xlim: xlim, defaults to [-2, 80]
    :type xlim: list, optional
    :param show: show figure
    :type show: bool
    :return: PlotR object
    :rtype: PlotR
    """
    pr = PlotR(rfsta, enf=enf, xlim=xlim, key=key)
    pr.plot(out_path=out_path, outformat=outformat, show=show)
    return pr

if __name__ == '__main__':
    pass