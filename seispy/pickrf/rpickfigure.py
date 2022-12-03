import glob
import obspy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from seispy.pickrf.pickfigure import RFFigure, indexpags
from seispy.rfcorrect import RFStation
from os.path import join, basename


class RPickFigure(RFFigure):
    def __init__(self, rfpath, width=8, height=10, dpi=100, xlim=[-2, 85]):
        super().__init__(rfpath, width=width, height=height, dpi=dpi, xlim=xlim)
        self.enf = 8
        
    def init_canvas(self, order='baz'):
        return super().init_canvas(order=order)

    def init_figure(self, width, height, dpi):
        self.fig = Figure(figsize=(width, height), dpi=dpi) 
        self.axr = self.fig.add_axes([0.15, 0.05, 0.55, 0.84])
        self.axb = self.fig.add_axes([0.755, 0.05, 0.22, 0.84])
        self.axg = self.axb.twiny()
    
    def read_sac(self, dt=0.1, order='baz'):
        if not isinstance(order, str):
            raise TypeError('The order must be str type')
        elif not order in ['baz', 'dis']:
            raise ValueError('The order must be \'baz\' or \'dis\'')
        else:
            pass
        if len(glob.glob(join(self.rfpath, '*P_R.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*_R.sac'))  
            self.comp = 'R'
        elif len(glob.glob(join(self.rfpath, '*P_Q.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*P_Q.sac'))
            self.comp = 'Q'
        elif len(glob.glob(join(self.rfpath, '*S_L.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*S_L.sac'))
            self.comp = 'L'
            self.enf = 4
        elif len(glob.glob(join(self.rfpath, '*S_Z.sac'))) != 0:
            tmp_files = glob.glob(join(self.rfpath, '*S_Z.sac'))
            self.comp = 'Z'
            self.enf = 4
        else:
            self.log.RFlog.info('No such file in the format of evt_phase_component.sac')
            exit()
        self.log.RFlog.info('Reading PRFs from {}'.format(self.rfpath))
        self.filenames = [basename(sac_file).split('_')[0] for sac_file in sorted(tmp_files)]
        self.phases = [basename(sac_file).split('_')[1] for sac_file in sorted(tmp_files)]
        self.rrf = obspy.read(join(self.rfpath, '*_{}.sac'.format(self.comp))).sort(['starttime']).resample(1/dt, window=None)
        self.time_axis = self.rrf[0].times() + self.rrf[0].stats.sac.b
        if self.xlim[1] > self.time_axis[-1]:
            self.xlim[1] = self.time_axis[-1]
        self.evt_num = len(self.rrf)
        self.log.RFlog.info('A total of {} PRFs loaded'.format(self.evt_num))
        self.baz = np.array([tr.stats.sac.baz for tr in self.rrf])
        self.gcarc = np.array([tr.stats.sac.gcarc for tr in self.rrf])
        self.rayp = np.array([tr.stats.sac.user0 for tr in self.rrf])
        self._sort(order)
        self.axpages, self.rfidx = indexpags(self.evt_num, self.maxidx)
        self.staname = (self.rrf[0].stats.network+'.'+self.rrf[0].stats.station).strip('.')
        self.fig.suptitle("%s (Latitude: %5.2f\N{DEGREE SIGN}, Longitude: %5.2f\N{DEGREE SIGN})" % (self.staname, self.rrf[0].stats.sac.stla, self.rrf[0].stats.sac.stlo), fontsize=20)

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
        # self.gcarc = [self.rrf[i].stats.sac.gcarc for i in range(self.evt_num)]
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
        self.set_ax_baz_dis()
    
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
        self.rlines = [[]]*self.evt_num
        self.rwvfillpos = [[]]*self.evt_num
        self.rwvfillnag = [[]]*self.evt_num
    
    def plotwave(self):
        bound = np.zeros(self.time_axis.shape[0])
        for i in range(self.evt_num):
            r_amp_axis = self.rrf[i].data*self.enf+i+1
            self.rlines[i], = self.axr.plot(self.time_axis, r_amp_axis, color="black", linewidth=0.2)            
            self.rwvfillpos[i] = self.axr.fill_between(self.time_axis, r_amp_axis, bound+i+1, where=r_amp_axis >i+1, facecolor='red', alpha=0.3)
            self.rwvfillnag[i] = self.axr.fill_between(self.time_axis, r_amp_axis, bound+i+1, where=r_amp_axis <i+1, facecolor='blue', alpha=0.3)

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
        goodidx = np.where(self.goodrf == 1)[0]
        newrfs = obspy.Stream([self.rrf[idx] for idx in goodidx])
        stadata = RFStation.read_stream(newrfs, self.rayp[goodidx],
                                        self.baz[goodidx], prime_comp=self.comp)
        stadata.event = np.array([self.filenames[i] for i in goodidx])
        self.plotfig = stadata.plotr(enf=self.enf, xlim=self.xlim)