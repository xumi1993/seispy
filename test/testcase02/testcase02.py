import numpy as np
from os.path import dirname, join
from obspy.io.sac import SACTrace
from seispy.hk import hkstack, ci, print_result
import re


def readpara(path):
    par_file = join(path, 'raysum-params')
    with open(par_file) as f:
        content = f.read()
    consp = content.split('\n')
    npts = int(consp[3])
    dt = float(consp[5])
    shift = float(consp[11])
    return npts, dt, shift


def initpara():
    h = np.arange(40, 80, 0.1)
    kappa = np.arange(1.6, 1.9, 0.01)
    path = join(dirname(__file__), 'benchmark')
    npts, dt, shift = readpara(path)
    baz, rayp = np.loadtxt(join(path, 'sample.geom'), usecols=(0, 1), unpack=True)
    rayp *= 1000
    ev_num = baz.shape[0]
    seis = np.empty([ev_num, npts])
    for _i, b in enumerate(baz):
        sac = SACTrace.read(join(path, 'tr_{:d}_{:.3f}.r'.format(int(b), rayp[0])))
        seis[_i] = sac.data
    return h, kappa, baz, rayp, npts, dt, shift, ev_num, seis


def subtc1(seis, shift, dt, rayp, h, kappa):
    try:
        stack, _, allstack, _ = hkstack(seis, shift, dt, rayp, h, kappa, vp=6.54)
        print('[Sub TC1 passed]: All PRFs have been stacked')
        return stack, allstack
    except Exception as e:
        raise RuntimeError('{}'.format(e))


def subtc2(allstack, h, kappa, ev_num):
    try:
        besth, bestk, cvalue, maxhsig, maxksig = ci(allstack, h, kappa, ev_num)
        print('[Sub TC2 passed]:')
        print_result(besth, bestk, maxhsig, maxksig, print_comment=True)
        return besth, bestk, cvalue, maxhsig, maxksig
    except Exception as e:
        raise RuntimeError('{}'.format(e))


if __name__ == '__main__':
    h, kappa, baz, rayp, npts, dt, shift, ev_num, seis = initpara()
    stack, allstack = subtc1(seis, shift, dt, rayp, h, kappa)
    subtc2(allstack, h, kappa, ev_num)
