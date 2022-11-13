import numpy as np
import re
from obspy.io.sac.sactrace import SACTrace
import matplotlib.pyplot as plt
from os.path import dirname, join
from seispy.rfcorrect import RFStation
from seispy.hkpara import hkpara
from seispy.geo import srad2skm
import argparse
from seispy.utils import load_cyan_map


def transarray(array, axis=0):
    if not isinstance(array, np.ndarray):
        raise ValueError('array should be `numpy.ndarray`')
    if len(array.shape) != 1:
        raise ValueError('array should be 1-d array')
    if axis == 0:
        return array.reshape(-1, array.shape[0])
    elif axis == 1:
        return array.reshape(array.shape[0], -1)
    else:
        raise ValueError('axis should be 0 or 1')


def vslow(v, rayp):
    return np.sqrt(1/(v**2) - rayp**2)


def tps(depth, eta_p, eta_s):
    return np.dot(transarray(eta_s - eta_p, axis=1), transarray(depth, axis=0))


def tppps(depth, eta_p, eta_s):
    return np.dot(transarray(eta_s + eta_p, axis=1), transarray(depth, axis=0))


def tpsps(depth, eta_s):
    return np.dot(transarray(2 * eta_s, axis=1), transarray(depth, axis=0))


def time2idx(times, ti0, dt):
    ti = ti0 + np.around(times / dt)
    return ti.reshape(ti.size).astype(int)


def hkstack(seis, t0, dt, p, h, kappa, vp=6.3, weight=(0.7, 0.2, 0.1)):
    # get dimensions
    nh = len(h)
    nk = len(kappa)
    nrf = len(p)

    # check the orientation of the seis array
    if seis.shape[0] != nrf:
        seis = seis.T
        if seis.shape[0] != nrf:
            raise IndexError('SEIS array dimensions should be (nt x nrf)')

    # amp correction for Ps
    am_cor = 151.5478 * p ** 2 + 3.2896 * p + 0.2618

    # get all vs, single column
    vs = vp / kappa

    # get index of direct P
    ti0 = round(t0 / dt)

    # initialize stacks
    tstack = np.zeros((nk, nh, 3))
    stack = np.zeros((nk, nh, 3))
    stack2 = np.zeros((nk, nh, 3))

    allstack = np.zeros((nk, nh, nrf))

    for i in range(nrf):
        eta_p = vslow(vp, p[i])
        eta_s = vslow(vs, p[i])

        # get times of Ps for all combinations of vs and H
        t1 = time2idx(tps(h, eta_p, eta_s), ti0, dt)
        t2 = time2idx(tppps(h, eta_p, eta_s), ti0, dt)
        t3 = time2idx(tpsps(h, eta_s), ti0, dt)

        tstack[:, :, 0] = am_cor[i] * seis[i, t1].reshape(nk, nh)
        tstack[:, :, 1] = am_cor[i] * seis[i, t2].reshape(nk, nh)
        tstack[:, :, 2] = -am_cor[i] * seis[i, t3].reshape(nk, nh)

        stack += tstack
        stack2 += tstack ** 2

        allstack[:, :, i] = weight[0] * tstack[:, :, 0] + weight[1] * tstack[:, :, 1] + weight[2] * tstack[:, :, 2]

    stack = stack / nrf
    stackvar = (stack2 - stack ** 2) / (nrf ** 2)

    allstackvar = np.var(allstack, axis=2)
    allstack = np.mean(allstack, axis=2)
    Normed_stack = allstack - np.min(allstack)
    Normed_stack = Normed_stack / np.max(Normed_stack)
    return stack, stackvar, Normed_stack, allstackvar


def plot(stack, allstack, h, kappa, besth, bestk, cvalue, cmap=load_cyan_map(), title=None, path=None):
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8), sharex='col', sharey='row')
    xlim = (h[0], h[-1])
    ylim = (kappa[0], kappa[-1])
    if title is not None:
        f.suptitle(title, fontsize='large')
    ax1.imshow(stack[:, :, 0], cmap=cmap, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], aspect='auto', origin='lower')
    ax1.set_ylabel('$V_P/V_S$')
    ax1.set_title('Ps')
    ax2.imshow(stack[:, :, 1], cmap=cmap, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], aspect='auto', origin='lower')
    ax2.set_title('PpPs')
    ax3.imshow(stack[:, :, 2], cmap=cmap, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], aspect='auto', origin='lower')
    ax3.set_title('PsPs+PpSs')
    ax3.set_xlabel('Moho depth (km)')
    ax3.set_ylabel('$V_P/V_S$')
    im = ax4.imshow(allstack, cmap=cmap, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], aspect='auto', origin='lower')
    ax4.plot(besth, bestk, color='red', marker='s', markerfacecolor='none')
    ax4.contour(allstack, [cvalue, 1], colors='k', extent=[xlim[0], xlim[1], ylim[0], ylim[1]], origin='lower')
    ax4.plot(xlim, [bestk, bestk], color='red', linestyle='--', linewidth=0.6)
    ax4.plot([besth, besth], ylim, color='red', linestyle='--', linewidth=0.6)
    ax4.set_xlabel('Moho depth (km)')

    plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9)
    _, yy, _, ww = ax4.get_position().bounds
    cax = plt.axes([0.93, yy, 0.016, ww])
    plt.colorbar(im, cax=cax)
    if path is None:
        plt.show()
    else:
        f.savefig(path, format='png', dpi=400, bbox_inches='tight')


def ci(allstack, h, kappa, ev_num):
    """
    Search best H and kappa from stacked matrix.
    Calculate error for H and kappa
    :param allstack: stacked HK matrix
    :param h: 1-D array of H
    :param kappa: 1-D array of kappa
    :param ev_num: event number
    :return:
    """
    [i, j] = np.unravel_index(allstack.argmax(), allstack.shape)
    bestk = kappa[i]
    besth = h[j]

    cvalue = 1 - np.std(allstack.reshape(allstack.size)) / np.sqrt(ev_num)
    cs = plt.contour(h, kappa, allstack, [cvalue])
    cs_path = cs.collections[0].get_paths()[0].vertices
    maxhsig = (np.max(cs_path[:, 0]) - np.min(cs_path[:, 0])) / 2
    maxksig = (np.max(cs_path[:, 1]) - np.min(cs_path[:, 1])) / 2
    plt.close()
    return besth, bestk, cvalue, maxhsig, maxksig


def print_result(besth, bestk, maxhsig, maxksig, print_comment=True):
    header = 'H\tH_error\tk\tk_error\n'
    if print_comment:
        msg = '{}{:.1f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(header, besth, maxhsig, bestk, maxksig)
    else:
        msg = '{:.1f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(besth, maxhsig, bestk, maxksig)
    print(msg)


def hksta(hpara, isplot=False, isdisplay=False):
    stadata = RFStation(hpara.rfpath, only_r=True)
    stack, _, allstack, _ = hkstack(stadata.datar, stadata.shift, stadata.sampling, srad2skm(stadata.rayp),
                                    hpara.hrange, hpara.krange, vp=hpara.vp, weight=hpara.weight)
    besth, bestk, cvalue, maxhsig, maxksig = ci(allstack, hpara.hrange, hpara.krange, stadata.ev_num)
    with open(hpara.hklist, 'a') as f:
        f.write('{}\t{:.3f}\t{:.3f}\t{:.1f}\t{:.2f}\t{:.2f}\t{:.3f}\n'.format(stadata.staname, stadata.stla, stadata.stlo,
                                                                              besth, maxhsig, bestk, maxksig))
    title = '{}\nMoho depth = ${:.1f}\pm{:.2f}$ km\n$V_P/V_S$ = ${:.2f}\pm{:.3f}$'.format(stadata.staname, besth,
                                                                                     maxhsig, bestk, maxksig)
    if isdisplay:
        print_result(besth, bestk, maxhsig, maxksig, print_comment=True)
    if isplot:
        img_path = join(hpara.hkpath, stadata.staname+'_Hk.png')
        plot(stack, allstack, hpara.hrange, hpara.krange, besth, bestk, cvalue, title=title, path=img_path)
    else:
        plot(stack, allstack, hpara.hrange, hpara.krange, besth, bestk, cvalue, title=title)


def hk():
    parser = argparse.ArgumentParser(description="HK stacking for single station")
    parser.add_argument('cfg_file', type=str, help='Path to HK configure file')
    parser.add_argument('-v', help='Display results to standard output',
                        dest='isdisplay', action='store_true')
    arg = parser.parse_args()
    hpara = hkpara(arg.cfg_file)
    hksta(hpara, isplot=True, isdisplay=arg.isdisplay)


def hktest():
    h = np.arange(40, 80, 0.1)
    kappa = np.arange(1.6, 1.9, 0.01)
    path = '/Users/xumj/Researches/YN_crust/dip_study/syn_hk/A1'
    baz, rayp = np.loadtxt(join(path, 'sample.geom'), usecols=(0, 1), unpack=True)
    rayp *= 1000
    with open(join(path, 'sample.tr')) as f:
        content = f.read()
    head_list = re.findall(r'#.*\n\s+(\d*)\s+(\d*)\s+(\d*.\d*)\s+(\d*)\s+(\d*.\d*)\n', content)[0]
    ev_num = int(head_list[0])
    npts = int(head_list[1])
    dt = float(head_list[2])
    shift = float(head_list[4])
    seis = np.empty([ev_num, npts])
    for _i, b, r in zip(range(ev_num), baz, rayp):
        sac = SACTrace.read(join(path, 'tr_{:d}_{:.3f}.r'.format(int(b), r)))
        seis[_i] = sac.data

    stack, _, allstack, _ = hkstack(seis, shift, dt, rayp, h, kappa, vp=6.3)
    besth, bestk, cvalue, maxhsig, maxksig = ci(allstack, h, kappa, ev_num)
    print_result(besth, bestk, maxhsig, maxksig, print_comment=True)
    plot(stack, allstack, h, kappa, besth, bestk, cvalue)

def hk_sta_test():
    hpara = hkpara('/Users/xumj/Researches/YNRF/hk.cfg')
    hksta(hpara, isplot=True)

if __name__ == '__main__':
    hk_sta_test()
