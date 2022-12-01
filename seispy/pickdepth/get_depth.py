from seispy.ccp3d import CCP3D
from seispy.geo import extrema
from seispy.signal import smooth
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator, AutoLocator
from matplotlib.backend_bases import MouseButton


def search_peak(tr, cpara, depmin, depmax):
    idx_all = extrema(tr)
    idx_ex = np.where(tr[idx_all] > 0)[0]
    peak = tr[idx_all[idx_ex]]
    peak_dep = idx_all[idx_ex] * cpara.stack_val + cpara.stack_range[0]
    idx_moho = np.where((peak_dep>depmin) & (peak_dep<depmax))[0]
    try:
        idx_moho_max = idx_moho[np.nanargmax(peak[idx_moho])]
        dep_moho = peak_dep[idx_moho_max]
    except:
        dep_moho = np.nan
    return dep_moho

class GoodDepth():
    def __init__(self, stack_data_path, logger, width=7.2, height=11, dpi=100, smooth=10) -> None:
        self.stack_data_path = stack_data_path
        self.logger = logger
        self.ccp_data = CCP3D.read_stack_data(stack_data_path)
        self.good_depth = pd.DataFrame(columns=['depth', 'amp', 'count', 'ci_low', 'ci_high'])
        self.bin_idx = 0
        self.smooth = 10
        self.smooth_val = int(self.smooth/self.ccp_data.cpara.stack_val)
        self.create_fig(width=width, height=height, dpi=dpi)
        
    def create_fig(self, width=7.2, height=11, dpi=100):
        self.fig = plt.figure(figsize=(width, height), dpi=dpi, constrained_layout=True)
        self.gs = gridspec.GridSpec(10, 8, figure=self.fig)
        self.ax_ew = self.fig.add_subplot(self.gs[0:2, 0:8])
        self.ax_ns = self.fig.add_subplot(self.gs[2:, 0:2])
        self.ax_stack = self.fig.add_subplot(self.gs[2:, 2:6])
        self.ax_count = self.fig.add_subplot(self.gs[2:, 6:8])

    def get_dep(self, depmin, depmax):
        self.depmin = depmin
        self.depmax = depmax
        if depmin > self.ccp_data.cpara.stack_range[-1] or \
           depmax < self.ccp_data.cpara.stack_range[0] or \
           depmax < depmin:
           raise ValueError('Depth range is out of stacking range.')
        for i, bin_stack in enumerate(self.ccp_data.stack_data):
            tr = smooth(bin_stack['mu'], self.smooth_val, 'hanning')
            dep_moho = search_peak(tr, self.ccp_data.cpara, depmin, depmax)
            if np.isnan(dep_moho):
                stackamp = np.nan
                count = np.nan
                ci_low = np.nan
                ci_high = np.nan
            else:
                idx = int((dep_moho - self.ccp_data.cpara.stack_range[0])/self.ccp_data.cpara.stack_val)
                stackamp = bin_stack['mu'][idx]
                count = bin_stack['count'][idx]
                ci_low = bin_stack['ci'][idx, 0]
                ci_high = bin_stack['ci'][idx, 1]
            this_df = pd.DataFrame([[dep_moho, stackamp, count, ci_low, ci_high]],
                                   columns=['depth', 'amp', 'count', 'ci_low', 'ci_high'])
            self.good_depth = pd.concat([self.good_depth, this_df], ignore_index=True)
    
    def get_adjcent(self, idx, val=5):
        self.val = val
        self.adjcent_ns = np.linspace(idx[0]-val, idx[0]+val, 2*val+1, dtype=np.uint8)
        self.adjcent_ew = np.linspace(idx[1]-val, idx[1]+val, 2*val+1, dtype=np.uint8)
        self.ns_lat = np.zeros(self.adjcent_ns.size)
        self.ew_lon = np.zeros(self.adjcent_ew.size)
        self.stack_ew = np.zeros([2*val+1, self.ccp_data.cpara.stack_range.size])
        self.stack_ns = np.zeros([2*val+1, self.ccp_data.cpara.stack_range.size])
        for i, idx_lat in enumerate(self.adjcent_ns):
            try:
                bin_idx = self.ccp_data.bin_map[idx_lat, idx[1]]
                self.stack_ns[i] = self.ccp_data.stack_data[bin_idx]['mu']
                self.ns_lat[i] = self.ccp_data.bin_loca[bin_idx, 0]
            except:
                self.stack_ns[i] = np.nan
                self.ns_lat[i] = np.nan
        for i, idx_lon in enumerate(self.adjcent_ew):
            try:
                bin_idx = self.ccp_data.bin_map[idx[0], idx_lon]
                self.stack_ew[i] = self.ccp_data.stack_data[bin_idx]['mu']
                self.ew_lon[i] = self.ccp_data.bin_loca[bin_idx, 1]
            except:
                self.stack_ew[i] = np.nan
                self.ew_lon[i] = np.nan

    def plot_bin(self, **kwargs):
        idx = np.where(self.ccp_data.bin_map==self.bin_idx)
        idx = [idx[0][0], idx[1][0]]
        if self.depmin < 100:
            self.plot_min = 0
        else:
            self.plot_min = self.depmin*0.85
        self.plot_max = self.depmax*1.15
        # self.this_depth = self.good_depth.iloc[self.bin_idx]['depth']
        self.get_adjcent(idx)
        self.init_fig(**kwargs)
        self.logger.PickDepthlog.info('{}/{}: Depth = {} km'.format(
            self.bin_idx+1, self.ccp_data.bin_loca.shape[0],
            self.good_depth.iloc[self.bin_idx]['depth']
        ))

    def init_fig(self):
        self.plot_ew(self.ax_ew)
        self.plot_ns(self.ax_ns)
        self.get_prot_moho()
        self.plot_stack(self.ax_stack, self.ax_count)
    
    def plot_ew(self, ax):
        ax.cla()
        ax.grid(color='gray', linestyle='--', linewidth=0.4, axis='y')
        for i, pos_lon in enumerate(self.ew_lon):
            # pos_lon = self.bin_mat[self.bin_loc_idx[0], idx_lon, 1]
            amp = self.stack_ew[i] * self.ccp_data.cpara.slide_val/50 + pos_lon
            ax.plot(amp, self.ccp_data.cpara.stack_range, linewidth=0.2, color='black')
            if not np.isnan(amp).all():
                if i == self.val:
                    ax.fill_betweenx(self.ccp_data.cpara.stack_range, amp, pos_lon, where=amp > pos_lon, facecolor='k', alpha=0.7)
                else:
                    ax.fill_betweenx(self.ccp_data.cpara.stack_range, amp, pos_lon, where=amp > pos_lon, facecolor='r', alpha=0.7)
        ax.set_ylim(self.plot_min, self.plot_max)
        # ax.set_yticks(np.arange(0, 100, 20))
        ax.set_ylabel('Depth (km)')
        ax.set_xlim(self.ccp_data.bin_loca[self.bin_idx, 1]-self.ccp_data.cpara.center_bin[-1]*self.val-0.1,
                    self.ccp_data.bin_loca[self.bin_idx, 1]+self.ccp_data.cpara.center_bin[-1]*self.val+0.1)
        ax.set_xlabel('Longitude ($^\circ$)')
        ax.set_title('Index: {}, Lat: {:.2f}$^\circ$, Lon: {:.2f}$^\circ$'.format(
            self.bin_idx+1, self.ccp_data.bin_loca[self.bin_idx][0], self.ccp_data.bin_loca[self.bin_idx][1]))
        ax.invert_yaxis()

    def plot_ns(self, ax):
        ax.cla()
        ax.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
        for i, pos_lat in enumerate(self.ns_lat):
            # pos_lat = self.bin_mat[idx_lat, self.bin_loc_idx[1], 0]
            amp = self.stack_ns[i] * self.ccp_data.cpara.slide_val/50 + pos_lat
            ax.plot(self.ccp_data.cpara.stack_range, amp, linewidth=0.2, color='black')
            if not np.isnan(amp).all():
                if i == self.val:
                    ax.fill_between(self.ccp_data.cpara.stack_range, amp, pos_lat, where=amp > pos_lat, facecolor='k', alpha=0.7)
                else:
                    ax.fill_between(self.ccp_data.cpara.stack_range, amp, pos_lat, where=amp > pos_lat, facecolor='red', alpha=0.7)
        ax.set_xlim(self.plot_min, self.plot_max)
        # ax.set_xticks(np.arange(0, 100, 20))
        ax.set_xlabel('Depth (km)')
        ax.set_ylim(self.ccp_data.bin_loca[self.bin_idx, 0]-self.ccp_data.cpara.center_bin[-1]*self.val-0.1,
                    self.ccp_data.bin_loca[self.bin_idx, 0]+self.ccp_data.cpara.center_bin[-1]*self.val+0.1)
        # ax.set_ylim(self.ns_lat[0]-0.1, self.ns_lat[-1]+0.1)
        ax.set_ylabel('Latitude ($^\circ$)')

    def get_prot_moho(self):
        # bin_stack = self.ccp.stack_data[i]
        self.mu = self.ccp_data.stack_data[self.bin_idx]['mu']
        self.ci = self.ccp_data.stack_data[self.bin_idx]['ci']
        if self.smooth is not None:
            self.mu = smooth(self.mu, self.smooth_val, 'hanning')
            self.ci[:, 0] = smooth(self.ci[:, 0], self.smooth_val, 'hanning')
            self.ci[:, 1] = smooth(self.ci[:, 1], self.smooth_val, 'hanning')
        if np.isnan(self.mu).all():
            self.prot_moho = [np.nan]
            self.maxpeak = [np.nan]
        else:
            self.ex_idx = extrema(self.mu)
            self.prot_moho = self.ccp_data.cpara.stack_range[self.ex_idx]
            #  self.prot_moho = self.prot_moho[np.where(self.mu[self.prot_moho] > 0)]
            # self.prot_moho.sort()
            self.maxpeak = np.nanmax(self.mu) + 0.1

    def plot_stack(self, ax, ax_c):
        ax.cla()
        ax_c.cla()
        maxamp = np.max(np.abs(self.mu))*1.5
        ax.plot([0, 0], [0, 100], color='k', linewidth=0.5)
        ax.plot(self.mu, self.ccp_data.cpara.stack_range)
        ax.plot(self.ci[:, 0], self.ccp_data.cpara.stack_range, lw=0.5, ls='--', color='black')
        ax.plot(self.ci[:, 1], self.ccp_data.cpara.stack_range, lw=0.5, ls='--', color='black')
        for extr in self.prot_moho:
            ax.plot([-maxamp, maxamp], [extr, extr], linewidth=2, color='tomato')
        self.pltdepth, = ax.plot([-maxamp, maxamp], [self.good_depth.iloc[self.bin_idx]['depth']]*2, linewidth=2, color='g')
        # ax.set_xlim([-self.maxpeak, self.maxpeak])
        ax.set_xlabel('Amplitude')
        ax.set_ylim([self.plot_min, self.plot_max])
        # ax.set_yticks(np.arange(self.ccp_data.cpara.stack_range[0], self.ccp_data.cpara.stack_range[-1], 10))
        ax.set_ylabel('Depth (km)')
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.invert_yaxis()

        ax_c.barh(self.ccp_data.cpara.stack_range, self.ccp_data.stack_data[self.bin_idx]['count'], height=0.8, facecolor='tan')
        ax_c.set_ylim([self.plot_min, self.plot_max])
        # ax_c.set_yticks(ax.get_yticks())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax_c.yaxis.set_minor_locator(AutoMinorLocator())
        ax_c.set_xlabel('Count')
        ax_c.invert_yaxis()

    def _get_next_bin(self):
        self.bin_idx += 1
        while True:
            # mu = self.ccp_data.stack_data[self.bin_idx]['mu']
            if self.bin_idx >= self.ccp_data.bin_loca.shape[0]:
                self.bin_idx = self.ccp_data.bin_loca.shape[0]-1
                break
            if np.isnan(self.good_depth.iloc[self.bin_idx]['depth']):
                self.bin_idx += 1
            else:
                break
    
    def _get_previous_bin(self):
        self.bin_idx -= 1
        while True:
            # mu = self.ccp_data.stack_data[self.bin_idx]['mu']
            if self.bin_idx < 0:
                self.bin_idx = 0
                break
            if np.isnan(self.good_depth.iloc[self.bin_idx]['depth']):
                self.bin_idx -= 1
            else:
                break
    
    def page_up(self, **kwargs):
        self._get_previous_bin()
        self.plot_bin(**kwargs)
    
    def page_down(self, **kwargs):
        self._get_next_bin()
        self.plot_bin(**kwargs)

    def on_click(self, event):
        if event.inaxes == self.ax_stack:
            idx = np.argmin(np.abs(self.prot_moho - event.ydata))
            sub_value = np.min(np.abs(self.prot_moho - event.ydata))
            if event.button is MouseButton.LEFT:
                if sub_value < 10:
                    self.pltdepth.set_ydata([self.prot_moho[idx]]*2)
                    self.pltdepth.set_visible(True)
                    self.logger.PickDepthlog.info('Set the depth to {} km'.format(self.prot_moho[idx]))
                    self.good_depth.iloc[self.bin_idx]['depth'] = self.prot_moho[idx]
            elif event.button is MouseButton.RIGHT:
                self.pltdepth.set_visible(False)
                self.logger.PickDepthlog.info('Set the depth to NaN')
                self.good_depth.iloc[self.bin_idx]['depth'] = np.nan

    def write(self, fname):
        bin_df = pd.DataFrame(self.ccp_data.bin_loca, columns=['lat', 'lon'])
        full_pd = pd.concat([bin_df, self.good_depth['depth'],
                  self.good_depth['ci_low'], self.good_depth['ci_high'],
                  self.good_depth['count']], axis=1)
        full_pd.to_csv(fname, sep=' ', header=False,
                       index=False, na_rep='nan',
                       float_format='%.4f')