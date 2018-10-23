import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from seispy.rfcorrect import SACStation
from seispy.rf import CfgParser
import argparse
import numpy as np
from os.path import join
import sys


def init_figure():
    h = plt.figure(figsize=(8, 10))
    gs = GridSpec(1, 3)
    gs.update(wspace=0.25)
    axr = plt.subplot(gs[0, 0:-1])
    axr.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
    axb = plt.subplot(gs[0, -1])
    axb.grid(color='gray', linestyle='--', linewidth=0.4, axis='x')
    return h, axr, axb


def read_process_data(lst):
    stadata = SACStation(lst, only_r=True)
    idx = np.argsort(stadata.bazi)
    stadata.event = stadata.event[idx]
    stadata.bazi = stadata.bazi[idx]
    stadata.datar = stadata.datar[idx]
    time_axis = np.arange(stadata.RFlength) * stadata.sampling - stadata.shift
    return stadata, time_axis


def plot_waves(axr, axb, stadata, time_axis, enf=12):
    bound = np.zeros(stadata.RFlength)
    for i in range(stadata.ev_num):
        datar = stadata.datar[i] * enf + (i + 1)
        # axr.plot(time_axis, stadata.datar[i], linewidth=0.2, color='black')
        axr.fill_between(time_axis, datar, bound + i+1, where=datar > i+1, facecolor='red',
                         alpha=0.7)
        axr.fill_between(time_axis, datar, bound + i+1, where=datar < i+1, facecolor='blue',
                         alpha=0.7)
    axb.scatter(stadata.bazi, np.arange(stadata.ev_num) + 1, s=7)


def set_fig(axr,  axb, stadata, station, xmin=-2, xmax=80):
    y_range = np.arange(stadata.ev_num) + 1
    x_range = np.arange(0, xmax+2, 5)
    space = 2

    # set axr
    axr.set_xlim(xmin, xmax)
    axr.set_xticks(x_range)
    axr.set_xticklabels(x_range, fontsize=8)
    axr.set_ylim(0, stadata.ev_num + space)
    axr.set_yticks(y_range)
    axr.set_yticklabels(stadata.event, fontsize=5)
    axr.set_xlabel('Time after P (s)', fontsize=13)
    axr.set_ylabel('Event', fontsize=13)
    axr.add_line(Line2D([0, 0], axr.get_ylim(), color='black'))
    axr.set_title('R components ({})'.format(station), fontsize=16)

    # set axb
    axb.set_xlim(0, 360)
    axb.set_xticks(np.linspace(0, 360, 7))
    axb.set_xticklabels(np.linspace(0, 360, 7, dtype='i'), fontsize=8)
    axb.set_ylim(0, stadata.ev_num + space)
    axb.set_yticks(y_range)
    axb.set_yticklabels(y_range, fontsize=5)
    axb.set_xlabel(r'Back-azimuth ($^\circ$)', fontsize=13)


def plotr(station, cfg_file, enf=6):
    pa = CfgParser(cfg_file)
    # pa.rfpath = join(pa.rfpath, station)
    lst = join(pa.rfpath, station + 'finallist.dat')
    h, axr, axb = init_figure()
    stadata, time_axis = read_process_data(lst)
    plot_waves(axr, axb, stadata, time_axis, enf=enf)
    set_fig(axr, axb, stadata, station)
    h.savefig(join(pa.imagepath, station + '_RT_bazorder_{:.1f}.pdf'.format(stadata.f0[0])), format='pdf')
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot R&T receiver functions")
    parser.add_argument('-s', help='Station as folder name of RFs and list', dest='station', type=str)
    parser.add_argument('-e', help='Enlargement factor', dest='enf', type=int, default=6)
    parser.add_argument('cfg_file', type=str, help='Path to configure file')
    arg = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    plotr(arg.station, arg.cfg_file, enf=arg.enf)


if __name__ == '__main__':
    station = 'XE.ES01'
    cfg_file = '/Users/xumj/Researches/Tibet_MTZ/process/paraRF.cfg'
    plotr(station, cfg_file)