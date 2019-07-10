import glob
import os
import obspy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from os.path import join, basename, dirname
from seispy.setuplog import setuplog
from seispy.plotRT import init_figure, set_fig, plot_waves


def indexpags(evt_num, maxidx=20):
    axpages = int(np.floor(evt_num/maxidx)+1)
    rfidx = []
    for i in range(axpages-1):
        rfidx.append(range(maxidx*i, maxidx*(i+1)))
    rfidx.append(range(maxidx*(axpages-1), evt_num))
    return axpages, rfidx


class StaData():
    def __init__(self, filenames, rrf, trf, baz, goodrf):
        goodidx = np.where(goodrf==1)
        self.event = [filenames[i] for i in goodidx[0]]
        self.bazi = baz[goodidx]
        self.RFlength = len(rrf[0].data)
        self.ev_num = len(self.bazi)
        self.datar = np.empty([self.ev_num, self.RFlength])
        self.datat = np.empty([self.ev_num, self.RFlength])
        for i, idx in enumerate(goodidx[0]):
            self.datar[i, :] = rrf[idx].data
            self.datat[i, :] = trf[idx].data


class RFFigure(Figure):
    def __init__(self, rfpath, width=21, height=11, dpi=100):
        super(RFFigure, self).__init__()
        
        self.log = setuplog()
        self.rfpath = rfpath
        self.enf = 3
        self.ipage = 0
        self.maxidx = 20
        self.xlim = [-2, 30]

        self.init_figure(width=21, height=11, dpi=100)
        self.read_sac()
        self.set_page()
        self.init_variables()
        self.plotwave()
        self.plotbaz()

    def set_ylabels(self):
        self.axr.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2)
        self.axr.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2))
        ylabels = self.filenames[self.rfidx[self.ipage][0]::]
        ylabels.insert(0, '')
        self.axr.set_yticklabels(ylabels)
        self.axt.set_ylim(self.axr.get_ylim())
        self.axt.set_yticks(self.axr.get_yticks())
        self.axb.set_ylim(self.axr.get_ylim())
        self.axb.set_yticks(self.axr.get_yticks())
        self.azi_label = ['%5.2f' % self.baz[i] for i in self.rfidx[self.ipage]]
        self.azi_label.insert(0, "")
        self.axb.set_yticklabels(self.azi_label)

    def set_page(self):
        self.set_ylabels()
        self.axr.plot([0, 0], [0, self.axpages*self.maxidx], color="black")
        self.axt.plot([0, 0], [0, self.axpages*self.maxidx], color="black")
        self.axr.set_xlim(self.xlim[0], self.xlim[1])
        self.axr.set_xticks(np.arange(self.xlim[0], self.xlim[1]+1, 2))
        self.axt.set_xlim(self.xlim[0], self.xlim[1])
        self.axt.set_xticks(np.arange(self.xlim[0], self.xlim[1]+1, 2))
        self.axb.set_xlim(0, 360)
        self.axb.set_xticks(np.arange(0, 361, 60))

    def init_variables(self):
        self.goodrf = np.ones(self.evt_num)
        self.rlines = [[] for i in range(self.evt_num)]
        self.tlines = [[] for i in range(self.evt_num)]
        self.rwvfillpos = [[] for i in range(self.evt_num)]
        self.rwvfillnag = [[] for i in range(self.evt_num)]
        self.twvfillpos = [[] for i in range(self.evt_num)]
        self.twvfillnag = [[] for i in range(self.evt_num)]

    def init_figure(self, width=21, height=11, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)  # 新建一个figure
        self.axr = self.fig.add_axes([0.1, 0.05, 0.35, 0.84])
        self.axt = self.fig.add_axes([0.47, 0.05, 0.35, 0.84])
        self.axb = self.fig.add_axes([0.855, 0.05, 0.12, 0.84])
        self.set_figure()

    def set_figure(self):
        self.axr.grid()
        self.axt.grid()
        self.axb.grid()
        self.axr.set_ylabel("Event")
        self.axr.set_xlabel("Time after P (s)")
        self.axr.set_title("R component")
        self.axt.set_xlabel("Time after P (s)")
        self.axt.set_title("T component")
        self.axb.set_xlabel("Backazimuth (\N{DEGREE SIGN})")

    def read_sac(self, dt=0.1):
        self.log.RFlog.info('Reading PRFs from {}'.format(self.rfpath))
        self.filenames = [basename(sac_file).split('_')[0] for sac_file in sorted(glob.glob(join(self.rfpath, '*_R.sac')))]
        self.rrf = obspy.read(join(self.rfpath, '*_R.sac')).sort(['starttime']).resample(1/dt)
        self.trf = obspy.read(join(self.rfpath, '*_T.sac')).sort(['starttime']).resample(1/dt)
        self.time_axis = self.rrf[0].times() + self.rrf[0].stats.sac.b
        self.evt_num = len(self.rrf)
        self.log.RFlog.info('A total of {} PRFs loaded'.format(self.evt_num))
        self.baz = np.array([tr.stats.sac.baz for tr in self.rrf])
        self.sort_baz_()
    
        self.axpages, self.rfidx = indexpags(self.evt_num, self.maxidx)
        self.staname = self.rrf[0].stats.network+'.'+self.rrf[0].stats.station
        self.fig.suptitle("%s (Latitude: %5.2f\N{DEGREE SIGN}, Longitude: %5.2f\N{DEGREE SIGN})" % (self.staname, self.rrf[0].stats.sac.stla, self.rrf[0].stats.sac.stlo), fontsize=20)

    def sort_baz_(self):
        idx = np.argsort(self.baz)
        self.baz = self.baz[idx]
        self.rrf = obspy.Stream([self.rrf[i] for i in idx])
        self.trf = obspy.Stream([self.trf[i] for i in idx])
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
        self.axb.scatter(self.baz, np.arange(self.evt_num)+1)

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

    def butprevious(self):
        self.ipage -= 1
        if self.ipage < 0:
            self.ipage = 0
            return
        self.set_ylabels()

    def butnext(self):
        self.ipage += 1
        if self.ipage >= self.axpages:
            self.ipage = self.axpages-1
            return
        if self.ipage == self.axpages-1:
            ymax = (self.ipage+1)*self.maxidx
            self.axr.set_ylim(self.rfidx[self.ipage][0], ymax)
            self.axr.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2))
            ylabels = self.filenames[self.rfidx[self.ipage][0]::]
            ylabels.insert(0, '')
            self.axr.set_yticklabels(ylabels)
            self.azi_label = ['%5.2f' % self.baz[i] for i in self.rfidx[self.ipage]]
        else:
            self.axr.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2)
            self.axr.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2))
            ylabels = self.filenames[self.rfidx[self.ipage][0]::]
            ylabels.insert(0, '')
            self.axr.set_yticklabels(ylabels)
            self.azi_label = ['%5.2f' % self.baz[i] for i in self.rfidx[self.ipage]]
        self.axt.set_ylim(self.axr.get_ylim())
        self.axt.set_yticks(self.axr.get_yticks())
        self.axb.set_ylim(self.axr.get_ylim())
        self.azi_label.insert(0, "")
        self.axb.set_yticklabels(self.azi_label)
        self.axb.set_yticks(self.axr.get_yticks())
    
    def enlarge(self):
        self.enf += 1
        self.axr.cla()
        self.axt.cla()
        self.plotwave()
        self.set_figure()
        self.set_page()

    def reduce(self):
        if self.enf > 1:
            self.enf -= 1
        else:
            self.enf = 1/(1/self.enf + 1)
        self.axr.cla()
        self.axt.cla()
        self.plotwave()
        self.set_figure()
        self.set_page()

    def finish(self):
        badidx = np.where(self.goodrf == 0)[0]
        self.log.RFlog.info("%d RFs are rejected" % len(badidx))
        with open(os.path.join(self.rfpath, basename(self.rfpath)+"finallist.dat"), 'w+') as fid:
            for i in range(self.evt_num):
                if self.goodrf[i] == 0:
                    files = glob.glob(join(self.rfpath, self.filenames[i]+'*.sac'))
                    for fil in files:
                        os.remove(fil)
                    self.log.RFlog.info("Reject PRF of "+self.filenames[i])
                    
                evla = self.rrf[i].stats.sac.evla
                evlo = self.rrf[i].stats.sac.evlo
                evdp = self.rrf[i].stats.sac.evdp
                dist = self.rrf[i].stats.sac.gcarc
                baz = self.rrf[i].stats.sac.baz
                rayp = self.rrf[i].stats.sac.user0
                mag = self.rrf[i].stats.sac.mag
                gauss = self.rrf[i].stats.sac.user1
                fid.write('%s %s %6.3f %6.3f %6.3f %6.3f %6.3f %8.7f %6.3f %6.3f\n' % (self.filenames[i], 'P', evla, evlo, evdp, dist, baz, rayp, mag, gauss))
    
    def plot(self, image_path):
        stadata = StaData(self.filenames, self.rrf, self.trf, self.baz, self.goodrf)
        h, axr, axt, axb, axr_sum, axt_sum = init_figure()
        plot_waves(axr, axt, axb, axr_sum, axt_sum, stadata, self.time_axis, enf=self.enf)
        set_fig(axr, axt, axb, axr_sum, axt_sum, stadata, self.staname)
        h.savefig(image_path)