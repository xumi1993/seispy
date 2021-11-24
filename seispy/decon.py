# -*- coding: utf-8 -*-
import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.util import next_pow_2
from math import pi
from scipy.fftpack import fft, ifft, ifftshift
from scipy.signal import fftconvolve, correlate
from scipy.linalg import solve_toeplitz
import matplotlib.pyplot as plt


def gaussFilter(dt, nft, f0):
    df = 1.0 / (nft * dt)
    nft21 = 0.5 * nft + 1
    f = df * np.arange(0, nft21)
    w = 2 * pi * f

    gauss = np.zeros([nft, 1])
    gauss1 = np.exp(-0.25 * (w / f0) ** 2) / dt
    gauss1.shape = (len(gauss1), 1)
    gauss[0:int(nft21)] = gauss1
    gauss[int(nft21):] = np.flipud(gauss[1:int(nft21) - 1])
    gauss = gauss[:, 0]

    return gauss


def gfilter(x, nfft, gauss, dt):
    Xf = fft(x, nfft)
    Xf = Xf * gauss * dt
    xnew = ifft(Xf, nfft).real
    return xnew


def correl(R, W, nfft):
    x = ifft(fft(R, nfft) * np.conj(fft(W, nfft)), nfft)
    x = x.real
    return x


def phaseshift(x, nfft, dt, tshift):
    Xf = fft(x, nfft)
    shift_i = int(tshift / dt)
    p = 2 * pi * np.arange(1, nfft + 1) * shift_i / nfft
    Xf = Xf * np.vectorize(complex)(np.cos(p), -np.sin(p))
    x = ifft(Xf, nfft) / np.cos(2 * pi * shift_i / nfft)
    x = x.real
    return x


def deconit(uin, win, dt, nt=None, tshift=10, f0=2.0, itmax=400, minderr=0.001, info=False, phase='P'):
    """
    Created on Wed Sep 10 14:21:38 2014
    [RFI, rms, it]=makeRFitdecon(uin,win,dt,nt,tshift,f0,itmax,minderr)

    In:
    uin = numerator (radial for PdS)
    win = denominator (vertical component for PdS)
    dt = sample interval (s)
    nt = number of samples
    tshift = Time until beginning of receiver function (s)
    f0 = width of gaussian filter
    itmax = max # iterations
    minderr = Min change in error required for stopping iterations

    Out:
    RFI = receiver function
    rms = Root mean square error for predicting numerator after each iteration

    @author: Mijian Xu @ NJU
    """
    # print('Iterative Decon (Ligorria & Ammon):\n')
    if len(uin) != len(win):
        raise ValueError('The two input trace must be in same length')
    elif nt is None:
        nt = len(uin)
    else:
        pass

    rms = np.zeros(itmax)
    nfft = next_pow_2(nt)
    p0 = np.zeros(nfft)

    u0 = np.zeros(nfft)
    w0 = np.zeros(nfft)

    u0[0:nt] = uin
    w0[0:nt] = win

    gaussF = gaussFilter(dt, nfft, f0)
    # gaussF = _gauss_filter(dt, nfft, f0)

    u_flt = gfilter(u0, nfft, gaussF, dt)
    w_flt = gfilter(w0, nfft, gaussF, dt)

    wf = fft(w0, nfft)
    r_flt = u_flt

    powerU = np.sum(u_flt ** 2)

    it = 0
    sumsq_i = 1
    d_error = 100 * powerU + minderr
    maxlag = 0.5 * nfft
    # print('\tMax Spike Display is ' + str((maxlag) * dt))

    while np.abs(d_error) > minderr and it < itmax:
        rw = correl(r_flt, w_flt, nfft)
        rw = rw / np.sum(w_flt ** 2)

        if phase == 'P':
            i1 = np.argmax(np.abs(rw[0:int(maxlag) - 1]))
        else:
            i1 = np.argmax(np.abs(rw))
        amp = rw[i1] / dt

        p0[i1] = p0[i1] + amp
        p_flt = gfilter(p0, nfft, gaussF, dt)
        p_flt = gfilter(p_flt, nfft, wf, dt)

        r_flt = u_flt - p_flt
        sumsq = np.sum(r_flt ** 2) / powerU
        rms[it] = sumsq
        d_error = 100 * (sumsq_i - sumsq)

        sumsq_i = sumsq

        it = it + 1

    p_flt = gfilter(p0, nfft, gaussF, dt)
    p_flt = phaseshift(p_flt, nfft, dt, tshift)
    RFI = p_flt[0:nt]
    rms = rms[0:it - 1]

    return RFI, rms, it 
    

def deconwater(uin, win, dt, tshift=10., wlevel=0.05, f0=2.0, normalize=False, phase='P'):
    """
    Frequency-domain deconvolution using waterlevel method.

    :param uin: R or Q component for the response function
    :type uin: np.ndarray
    :param win: Z or L component for the source function
    :type win: np.ndarray
    :param dt: sample interval in second 
    :type dt: float
    :param tshift: Time shift before P arrival, defaults to 10.
    :type tshift: float, optional
    :param wlevel: Waterlevel to stabilize the deconvolution, defaults to 0.05
    :type wlevel: float, optional
    :param f0: Gauss factor, defaults to 2.0
    :type f0: float, optional
    :param normalize: If normalize the amplitude of the RF, defaults to False
    :type normalize: bool, optional

    :return: (rf, rms) RF and final rms.
    :rtype: (np.ndarray, float)
    """
    if uin.size != win.size:
        raise ValueError('The length of the \'uin\' must be same as the \'win\'')
    nt = uin.size
    nft = next_pow_2(nt)
    nfpts = nft / 2 + 1     # number of freq samples
    fny = 1. / (2.* dt);     # nyquist
    delf = fny / (0.5 * nft)
    freq = delf * np.arange(nfpts)
    w = 2 * pi * freq

    # containers
    # rff = np.zeros(nft); # Rfn in freq domain
    upf = np.zeros(nft); # predicted numer in freq domain

    # Convert seismograms to freq domain
    uf = fft(uin, nft)
    wf = fft(win, nft)

    # denominator
    df = wf * wf.conjugate()
    dmax = max(df.real);   

    # add water level correction
    phi1 = wlevel * dmax # water level
    # nwl = length( find(df<phi1) ) # number corrected
    df[np.where(df.real < phi1)[0]] = phi1
    gaussF = gaussFilter(dt, nft, f0)
    nf = gaussF * uf * wf.conjugate()

    # compute RF
    rff = nf / df
    
    # compute predicted numerator
    upf = rff * wf

    # add phase shift to RF
    w = np.append(w, -np.flipud(w[1:-1]))
    rff = rff * np.exp(-1j * w * tshift)

    # back to time domain
    rft = ifft(rff , nft)
    rft = rft[0:nt].real

    # compute the fit
    uf = gaussF * uf # compare to filtered numerator
    ut = ifft(uf, nft).real
    upt = ifft(upf, nft).real

    powerU = np.sum(ut[0:nt] ** 2)
    rms = np.sum((upt[0:nt] - ut[0:nt]) ** 2 )/powerU

    if normalize:
        gnorm = np.sum(gaussF) * delf * dt
        rft = rft.real / gnorm

    return rft, rms, np.nan

def _add_zeros(a, numl, numr):
    """Add zeros at left and rigth side of array a"""
    return np.hstack([np.zeros(numl), a, np.zeros(numr)])


def _acorrt(a, num):
    """
    Not normalized auto-correlation of signal a.
    Sample 0 corresponds to zero lag time. Auto-correlation will consist of
    num samples. Correlation is performed in time domain with scipy.
    :param a: Data
    :param num: Number of returned data points
    :return: autocorrelation
    """
    return correlate(_add_zeros(a, 0, num - 1), a, 'valid')


def _xcorrt(a, b, num, zero_sample=0):
    """
    Not normalized cross-correlation of signals a and b.
    :param a,b: data
    :param num: The cross-correlation will consist of num samples.\n
        The sample with 0 lag time will be in the middle.
    :param zero_sample: Signals a and b are aligned around the middle of their
        signals.\n
        If zero_sample != 0 a will be shifted additionally to the left.
    :return: cross-correlation
    """
    if zero_sample > 0:
        a = _add_zeros(a, 2 * abs(zero_sample), 0)
    elif zero_sample < 0:
        a = _add_zeros(a, 0, 2 * abs(zero_sample))
    dif = len(a) - len(b) + 1 - num
    if dif > 0:
        b = _add_zeros(b, (dif + 1) // 2, dif // 2)
    else:
        a = _add_zeros(a, (-dif + 1) // 2, (-dif) // 2)
    return correlate(a, b, 'valid')


def decontime(uin, win, dt, tshift=10, spiking=1., normalize=True):
    nt = uin.size
    nshift = int(tshift * dt)
    STS = _acorrt(win, nt)
    STS = STS / STS[0]
    STS[0] += spiking
    STR = _xcorrt(uin, win, nt, nshift)
    assert len(STR) == len(STS)
    rf = solve_toeplitz(STS, STR)
    if normalize:
        norm = 1 / np.max(np.abs(rf))
        rf *= norm
    return rf


def deconvolute(uin, win, dt, method='iter', **kwargs):
    if method.lower() == 'iter':
        return deconit(uin, win, dt, **kwargs)
    elif method.lower() == 'water':
        return deconwater(uin, win, dt, **kwargs)
    else:
        raise ValueError('method must be in the \'iter\' or \'water\'')


if __name__ == '__main__':
    ldata = SACTrace.read('data/syn_S.L')
    qdata = SACTrace.read('data/syn_S.Q')
    l = ldata.data
    q = qdata.data
    # l = np.flip(ldata.data, axis=0)
    # q = np.flip(qdata.data, axis=0)
    time_axis = np.linspace(qdata.b, qdata.npts*qdata.delta+qdata.b, qdata.npts)

    rf, rms, it = deconit(l, q, qdata.delta, tshift=-qdata.b, f0=2, itmax=20, phase='S')
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10,6))
    axes[0].plot(time_axis, qdata.data)
    axes[1].plot(time_axis, ldata.data)
    axes[2].plot(time_axis, rf)
    plt.show()
