__author__ = 'xumj'


import matplotlib.pyplot as plt
import os, sys, glob
import obspy
import numpy as np
from matplotlib.widgets import Button
from operator import itemgetter
import initopts

def get_pos():
    (screen_width, screen_height) = plt.get_current_fig_manager().canvas.get_width_height()
    return screen_width, screen_height

def get_sac():
    path = '../53047'
    filesnames = glob.glob(path+'/*_R.sac')
    rffiles = obspy.read(path+'/*_R.sac')
    return rffiles.sort(['starttime']), filesnames

def indexpags(maxidx, evt_num):
    axpages = int(np.floor(evt_num/maxidx)+1)
    print(evt_num)
    rfidx = []
    for i in range(axpages-1):
        rfidx.append(range(maxidx*i, maxidx*(i+1)))
    rfidx.append(range(maxidx*(axpages-1), evt_num))
    return axpages, rfidx

class plotrffig():

    fig = plt.figure(figsize=(12, 12))
    axnext = plt.axes([0.81, 0.92, 0.07, 0.03])
    axprevious = plt.axes([0.71, 0.92, 0.07, 0.03])
    axfinish = plt.axes([0.91, 0.92, 0.07, 0.03])
    ax = plt.axes([0.1, 0.05, 0.6, 0.85])
    ax_baz = plt.axes([0.75, 0.05, 0.2, 0.85])
    ax.grid()
    ax_baz.grid()
    ax.set_ylabel("Event")
    ax.set_xlabel("Time after P (s)")
    ax_baz.set_xlabel("Backazimuth (\N{DEGREE SIGN})")
    bnext = Button(axnext, 'Next')
    bprevious = Button(axprevious, 'Previous')
    bfinish = Button(axfinish, 'Finish')

    def __init__(self, opts):
        self.opts = opts
        ax = self.ax
        ax_baz = self.ax_baz
        axpages, rfidx = indexpags(opts.maxidx, opts.evt_num)
        self.axpages = axpages
        self.rfidx = rfidx
        self.ipage = 0
        self.goodrf = np.ones(opts.evt_num)
        self.lines = [[] for i in range(opts.evt_num)]
        self.wvfillpos = [[] for i in range(opts.evt_num)]
        self.wvfillnag = [[] for i in range(opts.evt_num)]
        self.plotwave()
        self.plotbaz()
        ax.set_ylim(rfidx[self.ipage][0], rfidx[self.ipage][-1]+2)
        ax.set_yticks(np.arange(rfidx[self.ipage][0], rfidx[self.ipage][-1]+2))
        ax_baz.set_ylim(ax.get_ylim())
        ax_baz.set_yticks(ax.get_yticks())
        self.azi_label = ['%5.2f' % opts.baz[i] for i in rfidx[self.ipage]]
        self.azi_label.insert(0, "")
        ax_baz.set_yticklabels(self.azi_label)
        self.bnext.on_clicked(self.butnext)
        self.bprevious.on_clicked(self.butprevious)
        self.bfinish.on_clicked(self.finish)
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        ax.plot([0, 0], [0, opts.evt_num], color="black")
        ax.set_xlim(opts.xlim[0], opts.xlim[1])
        ax.set_xticks(np.arange(opts.xlim[0], opts.xlim[1]+1, 2))
        ax_baz.set_xlim(0, 360)
        ax_baz.set_xticks(np.arange(0, 361, 60))

    def finish(self, event):
        opts = self.opts
        badidx = np.where(self.goodrf == 0)[0]
        for i in badidx:
            print(opts.filenames[int(i)])
        sys.exit(0)

    def onclick(self, event):
        if not self.ax == event.inaxes:
            return
        click_idx = int(np.round(event.ydata))
        if self.goodrf[click_idx-1] == 1:
            self.goodrf[click_idx-1] = 0
            self.wvfillpos[click_idx-1].set_facecolor('gray')
            self.wvfillnag[click_idx-1].set_facecolor('gray')
        else:
            self.goodrf[click_idx-1] = 1
            self.wvfillpos[click_idx-1].set_facecolor('red')
            self.wvfillnag[click_idx-1].set_facecolor('blue')
        plt.draw()

    def plotwave(self):
        ax = self.ax
        opts = self.opts
        bound = np.zeros(opts.RFlength)
        time_axis = np.linspace(opts.b, opts.e, opts.RFlength)
        for i in opts.idx_bazi:
            rf = opts.rffiles[i]
            amp_axis = rf.data*opts.enf+i+1
            self.lines[i], = ax.plot(time_axis, amp_axis, color="black", linewidth=0.2)
            self.wvfillpos[i] = ax.fill_between(time_axis, amp_axis, bound+i+1, where=amp_axis >i+1, facecolor='red', alpha=0.3)
            self.wvfillnag[i] = ax.fill_between(time_axis, amp_axis, bound+i+1, where=amp_axis <i+1, facecolor='blue', alpha=0.3)

    def plotbaz(self):
        self.ax_baz.scatter(self.opts.baz, np.arange(self.opts.evt_num)+1)


    def butprevious(self, event):
        opts = self.opts
        ax = self.ax
        ax_baz = self.ax_baz
        self.ipage -= 1
        if self.ipage < 0:
            self.ipage = 0
            return
        ax.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2)
        ax.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2))
        ax_baz.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2)
        ax_baz.set_yticks(ax.get_yticks())
        self.azi_label = ['%5.2f' % opts.baz[i] for i in self.rfidx[self.ipage]]
        self.azi_label.insert(0, "")
        ax_baz.set_yticklabels(self.azi_label)
        plt.draw()

    def butnext(self, event):
        opts = self.opts
        ax = self.ax
        ax_baz = self.ax_baz
        self.ipage += 1
        if self.ipage >= self.axpages:
            self.ipage = self.axpages
            return
        ax.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2)
        ax.set_yticks(np.arange(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2))
        ax_baz.set_ylim(self.rfidx[self.ipage][0], self.rfidx[self.ipage][-1]+2)
        self.azi_label = ['%5.2f' % opts.baz[i] for i in self.rfidx[self.ipage]]
        self.azi_label.insert(0, "")
        ax_baz.set_yticklabels(self.azi_label)
        ax_baz.set_yticks(ax.get_yticks())
        plt.draw()

def main():
    opts = initopts.opts()
    opts.maxidx = 20
    opts.enf = 5
    opts.xlim = [-2, 30]
    opts.ylim = [0, 22]
    opts.rffiles, opts.filenames = get_sac()
    opts.evt_num = len(opts.rffiles)
    rf = opts.rffiles[0]
    dt = rf.stats.delta
    opts.b = rf.stats.sac.b
    opts.e = rf.stats.sac.e
    opts.RFlength = rf.data.shape[0]
    bazi = [tr.stats.sac.baz for tr in opts.rffiles]
    tmp_filenames = [[opts.filenames[i], bazi[i]] for i in range(opts.evt_num)]
    tmp_filenames = sorted(tmp_filenames, key=itemgetter(1))
    opts.filenames = [file[0] for file in tmp_filenames]
    opts.baz = np.sort(bazi)
    opts.idx_bazi = np.argsort(bazi)
    plotrf = plotrffig(opts)


if __name__ == "__main__":
    main()
    plt.show()
