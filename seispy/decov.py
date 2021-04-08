# -*- coding: utf-8 -*-
import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.util import next_pow_2
from math import pi
from scipy.fftpack import fft, ifft, ifftshift
from scipy.signal import fftconvolve
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


def decovit(uin, win, dt, nt=None, tshift=10, f0=2.0, itmax=400, minderr=0.001, info=False, phase='P'):
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

def __get_length(rsp_list):
    if isinstance(rsp_list, (list, tuple)):
        rsp_list = rsp_list[0]
    return len(rsp_list)


def _phase_shift_filter(nft, dt, tshift):
    """
    Construct filter to shift an array to account for time before onset

    :param nft: number of points for fft
    :param dt: sample spacing in seconds
    :param tshift: time to shift by in seconds
    :return: shifted array
    """
    freq = np.fft.fftfreq(nft, d=dt)
    return np.exp(-2j * pi * freq * tshift)


def _gauss_filter(dt, nft, f0, waterlevel=None):
    """
    Gaussian filter with width f0

    :param dt: sample spacing in seconds
    :param nft: length of filter in points
    :param f0: Standard deviation of the Gaussian Low-pass filter,
        corresponds to cut-off frequency in Hz for a response value of
        exp(0.5)=0.607.
    :param waterlevel: waterlevel for eliminating very low values
        (default: no waterlevel)
    :return: array with Gaussian filter frequency response
    """
    f = np.fft.fftfreq(nft, dt)
    gauss_arg = -0.5 * (f/f0) ** 2
    if waterlevel is not None:
        gauss_arg = np.maximum(gauss_arg, waterlevel)
    return np.exp(gauss_arg)


def deconv_waterlevel(rsp_list, src, sampling_rate, waterlevel=0.05, gauss=0.5,
                      tshift=10., length=None, normalize=0, nfft=None,
                      return_info=False):
    """
    Frequency-domain deconvolution using waterlevel method.

    Deconvolve src from arrays in rsp_list.

    :param rsp_list: either a list of arrays containing the response functions
        or a single array
    :param src: array with source function
    :param sampling_rate: sampling rate of the data
    :param waterlevel: waterlevel to stabilize the deconvolution
    :param gauss: Gauss parameter (standard deviation) of the
        Gaussian Low-pass filter,
        corresponds to cut-off frequency in Hz for a response value of
        exp(0.5)=0.607.
    :param tshift: delay time 0s will be at time tshift afterwards
    :param nfft: explicitely set number of samples for the fft
    :param length: number of data points in results, optional
    :param normalize: normalize all results so that the maximum of the trace
        with supplied index is 1. Set normalize to 'src' to normalize
        for the maximum of the prepared source. Set normalize to None for no
        normalization.
    :param return_info: return additionally a lot of different parameters in a
        dict for debugging purposes

    :return: (list of) array(s) with deconvolution(s)
    """
    if length is None:
        length = __get_length(rsp_list)
    N = length
    if nfft is None:
        nfft = next_pow_2(N)
    dt = 1. / sampling_rate
    ffilt = _phase_shift_filter(nfft, dt, tshift)
    if gauss is not None:
        ffilt = _gauss_filter(dt, nfft, gauss, waterlevel=-700) * ffilt
    spec_src = fft(src, nfft)
    spec_src_conj = np.conjugate(spec_src)
    spec_src_water = np.abs(spec_src * spec_src_conj)
    spec_src_water = np.maximum(
        spec_src_water, max(spec_src_water) * waterlevel)

    if normalize == 'src':
        spec_src = ffilt * spec_src * spec_src_conj / spec_src_water
        rf_src = ifft(spec_src, nfft)[:N]
        norm = 1 / max(rf_src)
        rf_src = norm * rf_src

    flag = False
    if not isinstance(rsp_list, (list, tuple)):
        flag = True
        rsp_list = [rsp_list]
    rf_list = [ifft(ffilt * fft(rsp, nfft) * spec_src_conj / spec_src_water,
                    nfft)[:N] for rsp in rsp_list]
    if normalize not in (None, 'src'):
        norm = 1. / max(rf_list[normalize])
    if normalize is not None:
        for rf in rf_list:
            rf *= norm
    if return_info:
        if normalize not in (None, 'src'):
            spec_src = ffilt * spec_src * spec_src_conj / spec_src_water
            rf_src = ifft(spec_src, nfft)[:N]
            norm = 1 / max(rf_src)
            rf_src = norm * rf_src
        info = {'rf_src': rf_src, 'rf_src_conj': spec_src_conj,
                'spec_src_water': spec_src_water,
                'freq': np.fft.fftfreq(nfft, d=dt),
                'gauss': ffilt, 'norm': norm, 'N': N, 'nfft': nfft}
        return rf_list, info
    elif flag:
        return rf
    else:
        return rf_list


if __name__ == '__main__':
    ldata = SACTrace.read('data/syn_S.L')
    qdata = SACTrace.read('data/syn_S.Q')
    l = ldata.data
    q = qdata.data
    # l = np.flip(ldata.data, axis=0)
    # q = np.flip(qdata.data, axis=0)
    time_axis = np.linspace(qdata.b, qdata.npts*qdata.delta+qdata.b, qdata.npts)

    rf, rms, it = decovit(l, q, qdata.delta, tshift=-qdata.b, f0=2, itmax=20, phase='S')
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10,6))
    axes[0].plot(time_axis, qdata.data)
    axes[1].plot(time_axis, ldata.data)
    axes[2].plot(time_axis, rf)
    plt.show()
