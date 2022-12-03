from curses import window
import glob
import os
import sys
import obspy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from os.path import join, basename
from seispy.setuplog import setuplog
from seispy.plotRT import init_figure, set_fig, plot_waves


def indexpags(evt_num, maxidx=20):
    full_pages = evt_num//maxidx
    if np.mod(evt_num, maxidx) == 0:
        axpages = full_pages
    else:
        axpages = full_pages+1
    rfidx = []
    for i in range(axpages-1):
        rfidx.append(np.arange(maxidx*i, maxidx*(i+1)))
    rfidx.append(np.arange(maxidx*(axpages-1), evt_num))
    return axpages, rfidx


class StaData():
    def __init__(self, filenames, rrf, trf, baz, goodrf):
        goodidx = np.where(goodrf == 1)
        self.event = [filenames[i] for i in goodidx[0]]
        self.bazi = baz[goodidx]
        self.rflength = len(rrf[0].data)
        self.ev_num = len(self.bazi)
        self.datar = np.empty([self.ev_num, self.rflength])
        self.datat = np.empty([self.ev_num, self.rflength])
        for i, idx in enumerate(goodidx[0]):    
            self.datar[i, :] = rrf[idx].data
            self.datat[i, :] = trf[idx].data
        self.time_axis = np.array([])


class RFFigure(Figure):
    def __init__(self, rfpath, width=21, height=11, dpi=100, xlim=[-2, 30]):
        super(RFFigure, self).__init__()
        self.width = width
        self.height = height
        self.dpi = dpi
        self.log = setuplog()
        self.rfpath = rfpath
        self.enf = 3.5
        self.ipage = 0
        self.maxidx = 20
        self.xlim = xlim
        self.plotfig = None

    def init_canvas(self, order='baz'):
        self.init_figure(width=self.width, height=self.height, dpi=self.dpi)
        self.read_sac(order=order)
        self.set_figure()
        self.set_page()
        self.init_variables()
        self.plotwave()
        self.plotbaz()

    def set_ylabels(self):
        self.axr.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][0]+self.maxidx+1)
        self.axr.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][0]+self.maxidx+1))
        ylabels = np.array(self.filenames)[self.rfidx[self.ipage]]
        ticklabels = [''] * (self.maxidx+1)        
        ticklabels[1: len(ylabels)+1] = ylabels
        self.axr.set_yticklabels(ticklabels)
        self.axt.set_ylim(self.axr.get_ylim())
        self.axt.set_yticks(self.axr.get_yticks())
        self.axb.set_ylim(self.axr.get_ylim())
        self.axb.set_yticks(self.axr.get_yticks())
        self.azi_label = ['%5.2f' % self.baz[i] for i in self.rfidx[self.ipage]]
        ticklabels = [''] * (self.maxidx+1)
        ticklabels[1: len(self.azi_label)+1] = self.azi_label
        self.axb.set_yticklabels(ticklabels)

    def set_ax_baz_dis(self):
        self.axb.set_xlim(0, 360)
        self.axb.set_xticks(np.arange(0, 361, 60))
        self.axb.xaxis.label.set_color('#1f77b4')
        self.axb.tick_params(axis='x', labelcolor='#1f77b4')
        self.axg.set_xlim(30, 90)
        self.axg.set_xticks(np.arange(30, 91, 10))
        self.axg.xaxis.label.set_color('#ff7f0e')
        self.axg.tick_params(axis='x', labelcolor='#ff7f0e')

    def set_page(self):
        self.set_ylabels()
        self.axr.plot([0, 0], [0, self.axpages*self.maxidx+1], color="black")
        self.axt.plot([0, 0], [0, self.axpages*self.maxidx+1], color="black")
        self.axr.set_xlim(self.xlim[0], self.xlim[1])
        self.axr.set_xticks(np.arange(self.xlim[0], self.xlim[1]+1, 2))
        self.axt.set_xlim(self.xlim[0], self.xlim[1])
        self.axt.set_xticks(np.arange(self.xlim[0], self.xlim[1]+1, 2))
        self.set_ax_baz_dis()

    def init_variables(self):
        self.goodrf = np.ones(self.evt_num)
        self.rlines = [[] for i in range(self.evt_num)]
        self.tlines = [[] for i in range(self.evt_num)]
        self.rwvfillpos = [[] for i in range(self.evt_num)]
        self.rwvfillnag = [[] for i in range(self.evt_num)]
        self.twvfillpos = [[] for i in range(self.evt_num)]
        self.twvfillnag = [[] for i in range(self.evt_num)]

    def init_figure(self, width=21, height=11, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi) 
        self.axr = self.fig.add_axes([0.1, 0.05, 0.35, 0.84])
        self.axt = self.fig.add_axes([0.47, 0.05, 0.35, 0.84])
        self.axb = self.fig.add_axes([0.855, 0.05, 0.12, 0.84])
        self.axg = self.axb.twiny()

    def set_figure(self):
        self.axr.grid(b=True, which='major', axis='x')
        self.axt.grid(b=True, which='major', axis='x')
        self.axb.grid(b=True, which='major')
        self.axr.set_ylabel("Event")
        self.axr.set_xlabel("Time after P (s)")
        self.axr.set_title("{} component".format(self.comp))
        self.axt.set_xlabel("Time after P (s)")
        self.axt.set_title("T component")
        self.axb.set_xlabel("Backazimuth (\N{DEGREE SIGN})")
        self.axg.set_xlabel('Distance (\N{DEGREE SIGN})')

    def read_sac(self, dt=0.1, order='baz'):
        if not isinstance(order, str):
            raise TypeError('The order must be str type')
        elif not order in ['baz', 'dis']:
            raise ValueError('The order must be \'baz\' or \'dis\'')
        else:
            pass
        if len(glob.glob(join(self.rfpath, '*P_R.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*P_R.sac'))  
            self.comp = 'R'
        elif len(glob.glob(join(self.rfpath, '*P_Q.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*P_Q.sac'))
            self.comp = 'Q'
        elif len(glob.glob(join(self.rfpath, '*S_L.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*S_L.sac'))
            self.comp = 'L'
        elif len(glob.glob(join(self.rfpath, '*S_Z.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*S_Z.sac'))
            self.comp = 'Z' 
        else:
            self.log.RFlog.error('No valid RFs in {}'.format(self.rfpath))
            sys.exit(1)
        self.log.RFlog.info('Reading PRFs from {}'.format(self.rfpath))
        self.filenames = [basename(sac_file).split('_')[0] for sac_file in sorted(tmp_files)]
        self.phases = [basename(sac_file).split('_')[1] for sac_file in sorted(tmp_files)]
        self.rrf = obspy.read(join(self.rfpath, '*_{}.sac'.format(self.comp))).sort(['starttime']).resample(1/dt, window=None)
        self.trf = obspy.read(join(self.rfpath, '*_T.sac')).sort(['starttime']).resample(1/dt, window=None)
        self.time_axis = self.rrf[0].times() + self.rrf[0].stats.sac.b
        self.evt_num = len(self.rrf)
        self.log.RFlog.info('A total of {} PRFs loaded'.format(self.evt_num))
        self.baz = np.array([tr.stats.sac.baz for tr in self.rrf])
        self.gcarc = np.array([tr.stats.sac.gcarc for tr in self.rrf])
        self._sort(order)
        self.axpages, self.rfidx = indexpags(self.evt_num, self.maxidx)
        self.staname = (self.rrf[0].stats.network+'.'+self.rrf[0].stats.station).strip('.')
        self.fig.suptitle("{} (Latitude: {:.2f}\N{DEGREE SIGN}, Longitude: {:.2f}\N{DEGREE SIGN})".format(
                          self.staname, self.rrf[0].stats.sac.stla, self.rrf[0].stats.sac.stlo), fontsize=20)

    def _sort(self, order):
        if order == 'baz':
            idx = np.argsort(self.baz)
        elif order == 'dis':
            idx = np.argsort(self.gcarc)
        else:
            pass
        self.baz = self.baz[idx]
        self.gcarc = self.gcarc[idx]
        self.phases = [self.phases[i] for i in idx]
        self.rrf = obspy.Stream([self.rrf[i] for i in idx])
        self.trf = obspy.Stream([self.trf[i] for i in idx])
        # self.gcarc = [self.rrf[i].stats.sac.gcarc for i in range(self.evt_num)]
        self.filenames = [self.filenames[i] for i in idx]

    def plotwave(self):
        bound = np.zeros(self.time_axis.shape[0])
        for i in range(self.evt_num):
            r_amp_axis = self.rrf[i].data*self.enf+i+1
            t_amp_axis = self.trf[i].data*self.enf+i+1
            self.rlines[i], = self.axr.plot(self.time_axis, r_amp_axis, color="black", linewidth=0.2)            
            self.tlines[i], = self.axt.plot(self.time_axis, t_amp_axis, color="black", linewidth=0.2)
            self.rwvfillpos[i] = self.axr.fill_between(self.time_axis, r_amp_axis, bound+i+1, where=r_amp_axis >i+1, facecolor='red', alpha=0.3)
            self.rwvfillnag[i] = self.axr.fill_between(self.time_axis, r_amp_axis, bound+i+1, where=r_amp_axis <i+1, facecolor='blue', alpha=0.3)
            self.twvfillpos[i] = self.axt.fill_between(self.time_axis, t_amp_axis, bound+i+1, where=t_amp_axis >i+1, facecolor='red', alpha=0.3)
            self.twvfillnag[i] = self.axt.fill_between(self.time_axis, t_amp_axis, bound+i+1, where=t_amp_axis <i+1, facecolor='blue', alpha=0.3)

    def plotbaz(self):
        self.axb.scatter(self.baz, np.arange(self.evt_num)+1, color='#1f77b4')
        self.axg.scatter(self.gcarc, np.arange(self.evt_num)+1, color='#ff7f0e', alpha=0.6)

    def onclick(self, event):
        if event.inaxes != self.axr and event.inaxes != self.axt:
            return
        click_idx = int(np.round(event.ydata))
        if click_idx > self.evt_num:
            return
        if self.goodrf[click_idx-1] == 1:
            self.log.RFlog.info("Selected "+self.filenames[click_idx-1])
            self.goodrf[click_idx-1] = 0
            self.rwvfillpos[click_idx-1].set_facecolor('gray')
            self.rwvfillnag[click_idx-1].set_facecolor('gray')
            self.twvfillpos[click_idx-1].set_facecolor('gray')
            self.twvfillnag[click_idx-1].set_facecolor('gray')
        else:
            self.log.RFlog.info("Canceled "+self.filenames[click_idx-1])
            self.goodrf[click_idx-1] = 1
            self.rwvfillpos[click_idx-1].set_facecolor('red')
            self.rwvfillnag[click_idx-1].set_facecolor('blue')
            self.twvfillpos[click_idx-1].set_facecolor('red')
            self.twvfillnag[click_idx-1].set_facecolor('blue')

    def _set_gray(self):
        for i in np.where(self.goodrf == 0)[0]:
            self.rwvfillpos[i].set_facecolor('gray')
            self.rwvfillnag[i].set_facecolor('gray')
            self.twvfillpos[i].set_facecolor('gray')
            self.twvfillnag[i].set_facecolor('gray')

    def butprevious(self):
        self.ipage -= 1
        if self.ipage < 0:
            self.ipage = 0
            return
        self.set_ylabels()
        self.set_figure()

    def butnext(self):
        self.ipage += 1
        if self.ipage >= self.axpages:
            self.ipage = self.axpages-1
            return
        self.set_ylabels()
        self.set_figure()

    def enlarge(self):
        self.enf += 1
        self.axr.cla()
        self.axt.cla()
        self.plotwave()
        self._set_gray()
        self.set_page()
        self.set_figure()

    def reduce(self):
        if self.enf > 1:
            self.enf -= 1
        else:
            self.enf = 1/(1/self.enf + 1)
        self.axr.cla()
        self.axt.cla()
        self.plotwave()
        self._set_gray()
        self.set_page()
        self.set_figure()

    def finish(self):
        badidx = np.where(self.goodrf == 0)[0]
        self.log.RFlog.info("%d RFs are rejected" % len(badidx))
        with open(os.path.join(self.rfpath, self.staname+"finallist.dat"), 'w+') as fid:
            for i in range(self.evt_num):
                if self.goodrf[i] == 0:
                    files = glob.glob(join(self.rfpath, self.filenames[i]+'*.sac'))
                    for fil in files:
                        os.remove(fil)
                    self.log.RFlog.info("Reject PRF of "+self.filenames[i])
                else:
                    evla = self.rrf[i].stats.sac.evla
                    evlo = self.rrf[i].stats.sac.evlo
                    evdp = self.rrf[i].stats.sac.evdp
                    rayp = self.rrf[i].stats.sac.user0
                    mag = self.rrf[i].stats.sac.mag
                    gauss = self.rrf[i].stats.sac.user1
                    fid.write('%s %s %6.3f %6.3f %6.3f %6.3f %6.3f %8.7f %6.3f %6.3f\n' % (
                              self.filenames[i], self.phases[i], evla, evlo, evdp, self.gcarc[i], self.baz[i], rayp, mag, gauss))

    def plot(self):
        plt.ion()
        plt.rcParams['toolbar'] = 'None'
        stadata = StaData(self.filenames, self.rrf, self.trf, self.baz, self.goodrf)
        stadata.time_axis = self.time_axis
        self.plotfig, axr, axt, axb, axr_sum, axt_sum = init_figure()
        plot_waves(axr, axt, axb, axr_sum, axt_sum, stadata, enf=self.enf)
        set_fig(axr, axt, axb, axr_sum, axt_sum, stadata, self.staname)
