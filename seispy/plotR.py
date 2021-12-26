import re
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from seispy.rfcorrect import RFStation
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


def read_process_data(rfpath):
    stadata = RFStation(rfpath, only_r=True)
    idx = np.argsort(stadata.bazi)
    stadata.event = stadata.event[idx]
    stadata.bazi = stadata.bazi[idx]
    stadata.datar = stadata.datar[idx]
    return stadata


def plot_waves(axr, axb, stadata, enf=12):
    bound = np.zeros(stadata.rflength)
    for i in range(stadata.ev_num):
        datar = stadata.datar[i] * enf + (i + 1)
        # axr.plot(time_axis, stadata.datar[i], linewidth=0.2, color='black')
        axr.fill_between(stadata.time_axis, datar, bound + i+1, where=datar > i+1, facecolor='red',
                         alpha=0.7)
        axr.fill_between(stadata.time_axis, datar, bound + i+1, where=datar < i+1, facecolor='blue',
                         alpha=0.7)
    axb.scatter(stadata.bazi, np.arange(stadata.ev_num) + 1, s=7)


def set_fig(axr,  axb, stadata, xmin=-2, xmax=80):
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
    axr.set_title('R components ({})'.format(stadata.staname), fontsize=16)

    # set axb
    axb.set_xlim(0, 360)
    axb.set_xticks(np.linspace(0, 360, 7))
    axb.set_xticklabels(np.linspace(0, 360, 7, dtype='i'), fontsize=8)
    axb.set_ylim(0, stadata.ev_num + space)
    axb.set_yticks(y_range)
    axb.set_yticklabels(y_range, fontsize=5)
    axb.set_xlabel(r'Back-azimuth ($^\circ$)', fontsize=13)


def plotr(rfpath, outpath='./', xlim=[-2, 80], enf=6, format='pdf'):
    h, axr, axb = init_figure()
    stadata = read_process_data(rfpath)
    plot_waves(axr, axb, stadata, enf=enf)
    set_fig(axr, axb, stadata, xlim[0], xlim[1])
    h.savefig(join(outpath, stadata.staname + '_R_bazorder_{:.1f}.{}'.format(
              stadata.f0[0], format)), format=format, dpi=500)


def main():
    parser = argparse.ArgumentParser(description="Plot R&T receiver functions")
    parser.add_argument('rfpath', help='Path to PRFs with a \'finallist.dat\' in it', type=str)
    parser.add_argument('-e', help='Enlargement factor, defaults to 6', dest='enf', type=float, default=6, metavar='enf')
    parser.add_argument('-o', help='Output path without file name, defaults to current path', dest='output', default='./', type=str, metavar='outpath')
    parser.add_argument('-t', help='Specify figure format. f = \'.pdf\', g = \'.png\', defaults to \'g\'',
                    dest='format', default='g', type=str, metavar='f|g')
    parser.add_argument('-x', help='The max time scale in sec, defaults to 85s', default=85, type=float, metavar='max_time')

    arg = parser.parse_args()
    if arg.format not in ('f', 'g'):
        raise ValueError('Error: The format must be in \'f\' and \'g\'')
    elif arg.format == 'g':
        fmt = 'png'
    elif arg.format == 'f':
        fmt = 'pdf'
    plotr(arg.rfpath, arg.output, enf=arg.enf, xlim=[-2, arg.x], format=fmt)



if __name__ == '__main__':
    pass