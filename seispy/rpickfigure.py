import glob
import os
import obspy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from seispy.pickfigure import RFFigure, indexpags
from os.path import join, basename
from seispy.setuplog import setuplog
from seispy.plotR import init_figure, set_fig, plot_waves


class StaData():
    def __init__(self, filenames, rrf, baz, goodrf):
        goodidx = np.where(goodrf == 1)
        self.event = [filenames[i] for i in goodidx[0]]
        self.bazi = baz[goodidx]
        self.rflength = len(rrf[0].data)
        self.ev_num = len(self.bazi)
        self.datar = np.empty([self.ev_num, self.rflength])
        for i, idx in enumerate(goodidx[0]):
            self.datar[i, :] = rrf[idx].data
        self.time_axis = np.array([])
        self.staname = ''


class RPickFigure(RFFigure):
    def __init__(self, rfpath, width=8, height=10, dpi=100, xlim=[-2, 85]):
        super().__init__(rfpath, width=width, height=height, dpi=dpi, xlim=xlim)
        self.enf = 8
        
    def init_canvas(self):
        return super().init_canvas()

    def init_figure(self, width, height, dpi):
        self.fig = Figure(figsize=(width, height), dpi=dpi) 
        self.axr = self.fig.add_axes([0.15, 0.05, 0.55, 0.84])
        self.axb = self.fig.add_axes([0.755, 0.05, 0.22, 0.84])
        self.axg = self.axb.twiny()
    
    def read_sac(self, dt=0.1):
        if len(glob.glob(join(self.rfpath, '*_R.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*_R.sac'))  
            self.comp = 'R'
        elif len(glob.glob(join(self.rfpath, '*_L.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*_L.sac'))
            self.comp = 'L'
        else:
            tmp_files = glob.glob(join(self.rfpath, '*_Z.sac'))
            self.comp = 'Z' 
        self.log.RFlog.info('Reading PRFs from {}'.format(self.rfpath))
        self.filenames = [basename(sac_file).split('_')[0] for sac_file in sorted(tmp_files)]
        self.rrf = obspy.read(join(self.rfpath, '*_{}.sac'.format(self.comp))).sort(['starttime']).resample(1/dt)
        self.time_axis = self.rrf[0].times() + self.rrf[0].stats.sac.b
        self.evt_num = len(self.rrf)
        self.log.RFlog.info('A total of {} PRFs loaded'.format(self.evt_num))
        self.baz = np.array([tr.stats.sac.baz for tr in self.rrf])
        self.sort_baz_()
        self.axpages, self.rfidx = indexpags(self.evt_num, self.maxidx)
        self.staname = (self.rrf[0].stats.network+'.'+self.rrf[0].stats.station).strip('.')
        self.fig.suptitle("%s (Latitude: %5.2f\N{DEGREE SIGN}, Longitude: %5.2f\N{DEGREE SIGN})" % (self.staname, self.rrf[0].stats.sac.stla, self.rrf[0].stats.sac.stlo), fontsize=20)

    
    def sort_baz_(self):
        idx = np.argsort(self.baz)
        self.baz = self.baz[idx]
        self.rrf = [self.rrf[i] for i in idx]
        self.gcarc = [self.rrf[i].stats.sac.gcarc for i in range(self.evt_num)]
        self.filenames = [self.filenames[i] for i in idx]

    def set_figure(self):
        self.axr.grid(b=True, which='major', axis='x')
        self.axb.grid(b=True, which='major')
        self.axr.set_ylabel("Event")
        self.axr.set_xlabel("Time after P (s)")
        self.axr.set_title("{} component".format(self.comp))
        self.axb.set_xlabel("Backazimuth (\N{DEGREE SIGN})")
        self.axg.set_xlabel('Distance (\N{DEGREE SIGN})')

    def set_page(self):
        self.set_ylabels()
        self.axr.plot([0, 0], [0, self.axpages*self.maxidx+1], color="black")
        self.axr.set_xlim(self.xlim[0], self.xlim[1])
        self.axr.set_xticks(np.arange(self.xlim[0], self.xlim[1]+1, 2))
        self.axb.set_xlim(0, 360)
        self.axb.set_xticks(np.arange(0, 361, 60))
        self.axg.set_xlim(30, 90)
        self.axg.set_xticks(np.arange(30, 91, 10))
    
    def set_ylabels(self):
        self.axr.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][0]+self.maxidx+1)
        self.axr.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][0]+self.maxidx+1))
        ylabels = np.array(self.filenames)[self.rfidx[self.ipage]]
        ticklabels = [''] * (self.maxidx+1)        
        ticklabels[1: len(ylabels)+1] = ylabels
        self.axr.set_yticklabels(ticklabels)
        self.axb.set_ylim(self.axr.get_ylim())
        self.axb.set_yticks(self.axr.get_yticks())
        self.azi_label = ['%5.2f' % self.baz[i] for i in self.rfidx[self.ipage]]
        ticklabels = [''] * (self.maxidx+1)
        ticklabels[1: len(self.azi_label)+1] = self.azi_label
        self.axb.set_yticklabels(ticklabels)
    
    def init_variables(self):
        self.goodrf = np.ones(self.evt_num)
        self.rlines = [[] for i in range(self.evt_num)]
        self.rwvfillpos = [[] for i in range(self.evt_num)]
        self.rwvfillnag = [[] for i in range(self.evt_num)]
    
    def plotwave(self):
        bound = np.zeros(self.time_axis.shape[0])
        for i in range(self.evt_num):
            r_amp_axis = self.rrf[i].data*self.enf+i+1
            self.rlines[i], = self.axr.plot(self.time_axis, r_amp_axis, color="black", linewidth=0.2)            
            self.rwvfillpos[i] = self.axr.fill_between(self.time_axis, r_amp_axis, bound+i+1, where=r_amp_axis >i+1, facecolor='red', alpha=0.3)
            self.rwvfillnag[i] = self.axr.fill_between(self.time_axis, r_amp_axis, bound+i+1, where=r_amp_axis <i+1, facecolor='blue', alpha=0.3)

    def plotbaz(self):
        self.axb.scatter(self.baz, np.arange(self.evt_num)+1)
        self.axg.scatter(self.gcarc, np.arange(self.evt_num)+1, color='orange', alpha=0.6)

    def onclick(self, event):
        if event.inaxes != self.axr:
            return
        click_idx = int(np.round(event.ydata))
        if click_idx > self.evt_num:
            return
        if self.goodrf[click_idx-1] == 1:
            self.log.RFlog.info("Selected "+self.filenames[click_idx-1])
            self.goodrf[click_idx-1] = 0
            self.rwvfillpos[click_idx-1].set_facecolor('gray')
            self.rwvfillnag[click_idx-1].set_facecolor('gray')
        else:
            self.log.RFlog.info("Canceled "+self.filenames[click_idx-1])
            self.goodrf[click_idx-1] = 1
            self.rwvfillpos[click_idx-1].set_facecolor('red')
            self.rwvfillnag[click_idx-1].set_facecolor('blue')

    def _set_gray(self):
        for i in np.where(self.goodrf == 0)[0]:
            self.rwvfillpos[i].set_facecolor('gray')
            self.rwvfillnag[i].set_facecolor('gray')
    
    def enlarge(self):
        self.enf += 1
        self.axr.cla()
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
        self.plotwave()
        self._set_gray()
        self.set_page()
        self.set_figure()

    def plot(self):
        plt.ion()
        plt.rcParams['toolbar'] = 'None'
        stadata = StaData(self.filenames, self.rrf, self.baz, self.goodrf)
        stadata.time_axis = self.time_axis
        stadata.staname = self.staname
        self.plotfig, axr, axb = init_figure()
        plot_waves(axr, axb, stadata, enf=self.enf)
        set_fig(axr, axb, stadata, self.xlim[0], self.xlim[1])