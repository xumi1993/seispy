import numpy as np
from os.path import join
import re
from obspy.io.sac.sactrace import SACTrace


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
    return stack, stackvar, allstack, allstackvar


def hktest():
    h = np.arange(30, 70, 1)
    kappa = np.arange(1.6, 1.9, 0.01)
    path = '/Users/xumj/Researches/YN_crust/dip_study/syn_hk/benchmark'
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

    _, _, allstack, _ = hkstack(seis, shift, dt, rayp, h, kappa, vp=6.54)
    [i, j] = np.unravel_index(allstack.argmax(), allstack.shape)
    print('Moho depth = {:.0f}\nkappa = {:.2f}'.format(h[j], kappa[i]))


if __name__ == '__main__':
    hktest()


