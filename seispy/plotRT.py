import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import numpy as np
from os.path import join


class PlotRT:
    def __init__(self, stadata, enf=3, xlim=[2, 80], key='bazi') -> None:
        """Plot receiver function in R and T components

        :param stadata: receiver function data
        :type stadata: seispy.rfcorrect.RFStation
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
        self.fig = plt.figure(figsize=(11.7, 8.3))
        gs = GridSpec(17, 3)
        gs.update(wspace=0.25)
        self.axr_sum = plt.subplot(gs[0, 0])
        self.axr_sum.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        self.axr = plt.subplot(gs[1:, 0])
        self.axr.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        self.axt_sum = plt.subplot(gs[0, 1])
        self.axt_sum.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        self.axt = plt.subplot(gs[1:, 1])
        self.axt.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        self.axb = plt.subplot(gs[1:, 2])
        self.axb.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')

    def plot_waves(self):
        """Plot PRFs with R and T components
        """

        bound = np.zeros(self.stadata.rflength)
        for i in range(self.stadata.ev_num):
            datar = self.stadata.data_prime[i] * self.enf + (i + 1)
            datat = self.stadata.datat[i] * self.enf + (i + 1)
            # axr.plot(time_axis, stadata.datar[i], linewidth=0.2, color='black')
            self.axr.fill_between(self.stadata.time_axis, datar, bound + i+1, where=datar > i+1, facecolor='red',
                            alpha=0.7)
            self.axr.fill_between(self.stadata.time_axis, datar, bound + i+1, where=datar < i+1, facecolor='blue',
                            alpha=0.7)
            # axt.plot(time_axis, stadata.datat[i], linewidth=0.2, color='black')
            self.axt.fill_between(self.stadata.time_axis, datat, bound + i + 1, where=datat > i+1, facecolor='red',
                            alpha=0.7)
            self.axt.fill_between(self.stadata.time_axis, datat, bound + i + 1, where=datat < i+1, facecolor='blue',
                            alpha=0.7)
        datar = np.mean(self.stadata.data_prime, axis=0)
        datar /= np.max(datar)
        datat = np.mean(self.stadata.datat, axis=0)
        datat /= np.max(datar)
        self.axr_sum.fill_between(self.stadata.time_axis, datar, bound, where=datar > 0, facecolor='red', alpha=0.7)
        self.axr_sum.fill_between(self.stadata.time_axis, datar, bound, where=datar < 0, facecolor='blue', alpha=0.7)
        self.axt_sum.fill_between(self.stadata.time_axis, datat, bound, where=datat > 0, facecolor='red', alpha=0.7)
        self.axt_sum.fill_between(self.stadata.time_axis, datat, bound, where=datat < 0, facecolor='blue', alpha=0.7)
        self.axb.scatter(self.stadata.bazi, np.arange(self.stadata.ev_num) + 1, s=7)


    def set_fig(self):
        """Set figure
        """
        y_range = np.arange(self.stadata.ev_num) + 1
        x_range = np.arange(0, self.xlim[1]+2, 2)
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

        # set axr_sum
        self.axr_sum.set_title('{} components ({})'.format(self.stadata.comp, self.stadata.staname), fontsize=16)
        self.axr_sum.set_xlim(self.xlim)
        self.axr_sum.set_xticks(x_range)
        self.axr_sum.set_xticklabels([])
        self.axr_sum.set_ylim(-0.5, 1.25)
        self.axr_sum.set_yticks([0.375])
        self.axr_sum.set_yticklabels(['Sum'], fontsize=8)
        self.axr_sum.tick_params(axis='y', left=False)
        self.axr_sum.add_line(Line2D([0, 0], self.axr_sum.get_ylim(), color='black'))

        # set axt
        self.axt.set_xlim(self.xlim)
        self.axt.set_xticks(x_range)
        self.axt.set_xticklabels(x_range, fontsize=8)
        self.axt.set_ylim(0, self.stadata.ev_num + space)
        self.axt.set_yticks(y_range)
        bazi = ['{:.1f}'.format(ba) for ba in self.stadata.bazi]
        self.axt.set_yticklabels(bazi, fontsize=5)
        self.axt.set_xlabel('Time after P (s)', fontsize=13)
        self.axt.set_ylabel(r'Back-azimuth ($\circ$)', fontsize=13)
        self.axt.add_line(Line2D([0, 0], self.axt.get_ylim(), color='black'))

        # set axt_sum
        self.axt_sum.set_title('T components ({})'.format(self.stadata.staname), fontsize=16)
        self.axt_sum.set_xlim(self.xlim)
        self.axt_sum.set_xticks(x_range)
        self.axt_sum.set_xticklabels([])
        self.axt_sum.set_ylim(-0.5, 1.25)
        self.axt_sum.set_yticks([0.375])
        self.axt_sum.set_yticklabels(['Sum'], fontsize=8)
        self.axt_sum.tick_params(axis='y', left=False)
        self.axt_sum.add_line(Line2D([0, 0], self.axt_sum.get_ylim(), color='black'))

        # set axb
        self.axb.set_xlim(0, 360)
        self.axb.set_xticks(np.linspace(0, 360, 7))
        self.axb.set_xticklabels(np.linspace(0, 360, 7, dtype='i'), fontsize=8)
        self.axb.set_ylim(0, self.stadata.ev_num + space)
        self.axb.set_yticks(y_range)
        self.axb.set_yticklabels(y_range, fontsize=5)
        self.axb.set_xlabel(r'Back-azimuth ($\circ$)', fontsize=13)


    def plot(self, out_path=None, outformat='g', show=False):
        """Plot PRFs with R and T components

        :param rfsta: Path to PRFs
        :type rfsta: seispy.rfcorrect.RFStation
        :param enf: The enlarge factor, defaults to 3
        :type enf: int, optional
        :param out_path: The output path, defaults to current directory
        :type out_path: str, optional
        :param key: The key to sort PRFs, avialible for ``event``, ``evla``, ``evlo``, ``evdp``,
            ``dis``, ``bazi``, ``rayp``, ``mag``, ``f0``, defaults to ``bazi``
        :type key: str, optional
        :param outformat: File format of the image file, g as \'png\', f as \'pdf\', defaults to 'g'
        :type outformat: str, optional
        """
        self.plot_waves()
        self.set_fig()
        if out_path is not None:
            if outformat == 'g':
                self.fig.savefig(join(out_path, self.stadata.staname+'_{}T_{}order_{:.1f}.png'.format(
                                self.stadata.comp, self.key, self.stadata.f0[0])),
                                dpi=400, bbox_inches='tight')
            elif outformat == 'f':
                self.fig.savefig(join(out_path, self.stadata.staname+'_{}T_{}order_{:.1f}.pdf'.format(
                                self.stadata.comp, self.key, self.stadata.f0[0])),
                                format='pdf', bbox_inches='tight')
        if show:
            plt.show()

def plotrt(rfsta, out_path='./', xlim=[-2, 30], key='bazi', enf=6, outformat='g', show=False):
    """Plot PRFs with R and T components

    :param rfsta: Path to PRFs
    :type rfsta: seispy.rfcorrect.RFStation
    :param enf: The enlarge factor, defaults to 6
    :type enf: int, optional
    :param out_path: The output path, defaults to current directory
    :type out_path: str, optional
    :param key: The key to sort PRFs, avialible for ``event``, ``evla``, ``evlo``, ``evdp``,
        ``dis``, ``bazi``, ``rayp``, ``mag``, ``f0``, defaults to ``bazi``
    :type key: str, optional
    :param outformat: File format of the image file, g as \'png\', f as \'pdf\', defaults to 'g'
    :type outformat: str, optional
    :param xlim: xlim, defaults to [-2, 30]
    :type xlim: list, optional
    :param show: show figure
    :type show: bool
    """
    prt = PlotRT(rfsta, enf=enf, xlim=xlim, key=key)
    prt.plot(out_path, outformat, show)
    return prt


if __name__ == '__main__':
    pass
