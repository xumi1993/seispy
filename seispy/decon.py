# -*- coding: utf-8 -*-
import numpy as np
import obspy
from obspy.signal.util import next_pow_2
from numpy.fft import fft, ifft
# from scipy.linalg import solve_toeplitz


def gaussFilter(dt, nft, f0):
    """
    Gaussian filter in frequency domain.

    :param dt: sample interval in second
    :type dt: float
    :param nft: number of samples
    :type nft: int
    :param f0: Gauss factor
    :type f0: float

    :return: Gaussian filter in frequency domain
    :rtype: np.ndarray
    """
    df = 1.0 / (nft * dt)
    nft21 = 0.5 * nft + 1
    f = df * np.arange(0, nft21)
    w = 2 * np.pi * f

    gauss = np.zeros([nft, 1])
    gauss1 = np.exp(-0.25 * (w / f0) ** 2) / dt
    gauss1.shape = (len(gauss1), 1)
    gauss[0:int(nft21)] = gauss1
    gauss[int(nft21):] = np.flipud(gauss[1:int(nft21) - 1])
    gauss = gauss[:, 0]

    return gauss


def gfilter(x, nfft, gauss, dt):
    """
    Apply Gaussian filter on time series.

    :param x: input trace
    :type x: np.ndarray
    :param nfft: number of samples
    :type nfft: int
    :param gauss: Gaussian filter in frequency domain
    :type gauss: np.ndarray
    :param dt: sample interval in second
    :type dt: float

    :return: Filtered data in time domain
    :rtype: np.ndarray
    """
    Xf = fft(x, nfft)
    Xf = Xf * gauss * dt
    xnew = ifft(Xf, nfft).real
    return xnew


def correl(R, W, nfft):
    """
    Correlation in frequency domain.

    :param R: numerator
    :type R: np.ndarray
    :param W: denominator
    :type W: np.ndarray

    :return: Correlation in frequency domain
    :rtype: np.ndarray
    """
    x = ifft(fft(R, nfft) * np.conj(fft(W, nfft)), nfft)
    x = x.real
    return x


def phaseshift(x, nfft, dt, tshift):
    """
    Phase shift in frequency domain.
    
    :param x: input trace
    :type x: np.ndarray
    :param nfft: number of samples
    :type nfft: int
    :param dt: sample interval in second
    :type dt: float
    :param tshift: Time shift before P arrival
    :type tshift: float

    :return: Phase shifted data in time domain
    :rtype: np.ndarray
    """
    Xf = fft(x, nfft)
    shift_i = int(tshift / dt)
    p = 2 * np.pi * np.arange(1, nfft + 1) * shift_i / nfft
    Xf = Xf * np.vectorize(complex)(np.cos(p), -np.sin(p))
    x = ifft(Xf, nfft) / np.cos(2 * np.pi * shift_i / nfft)
    x = x.real
    return x


def deconit(uin, win, dt, nt=None, tshift=10, f0=2.0, itmax=400, minderr=0.001, phase='P'):
    """
    Iterative deconvolution using Ligorria & Ammon method.
    @author: Mijian Xu @ NJU
    Created on Wed Sep 10 14:21:38 2014 

    :param uin: R or Q component for the response function
    :type uin: np.ndarray
    :param win: Z or L component for the source function
    :type win: np.ndarray
    :param dt: sample interval in second
    :type dt: float
    :param nt: number of samples, defaults to None
    :type nt: int, optional
    :param tshift: Time shift before P arrival, defaults to 10.
    :type tshift: float, optional
    :param f0: Gauss factor, defaults to 2.0
    :type f0: float, optional
    :param itmax: Max iterations, defaults to 400
    :type itmax: int, optional
    :param minderr: Min change in error required for stopping iterations, defaults to 0.001
    :type minderr: float, optional
    :param phase: Phase of the RF, defaults to 'P'
    :type phase: str, optional

    :return: (RFI, rms, it) RF, rms and number of iterations.
    :rtype: (np.ndarray, np.ndarray, int)
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
    w = 2 * np.pi * freq

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

    return rft, rms


def deconvolute(uin, win, dt, method='iter', **kwargs):
    """ Deconvolute receiver function from waveforms.
    :param uin: R or Q component for the response function
    :type uin: np.ndarray
    :param win: Z or L component for the source function
    :type win: np.ndarray
    :param dt: sample interval in second
    :type dt: float
    :param method: Method for deconvolution, defaults to 'iter'
    :type method: str, optional
    :param kwargs: Parameters for deconvolution
    :type kwargs: dict
    :return: RF, rms, [iter]
    :rtype: np.ndarray, float, int
    """
    if method.lower() == 'iter':
        return deconit(uin, win, dt, **kwargs)
    elif method.lower() == 'water':
        return deconwater(uin, win, dt, **kwargs)
    else:
        raise ValueError('method must be \'iter\' or \'water\'')


class RFTrace(obspy.Trace):
    """ 
    Class for receiver function trace.
    """
    def __init__(self, data=..., header=None):
        super().__init__(data=data, header=header)

    @classmethod
    def deconvolute(cls, utr, wtr, method='iter', **kwargs):
        """
        Deconvolute receiver function from waveforms.

        :param utr: R or Q component for the response function
        :type utr: obspy.Trace
        :param wtr: Z or L component for the source function
        :type wtr: obspy.Trace
        :param method: Method for deconvolution, defaults to 'iter'
        :type method: str, optional
        :param kwargs: Parameters for deconvolution
        :type kwargs: dict

        :return: RFTrace object
        :rtype: RFTrace
        """
        header = utr.stats.__getstate__()
        for key, value in kwargs.items():
            header[key] = value
        if method.lower() == 'iter':
            rf, rms, it = deconit(utr.data, wtr.data, utr.stats.delta, **kwargs)
            header['rms'] = rms
            header['iter'] = it
        elif method.lower() == 'water':
            rf, rms = deconwater(utr.data, wtr.data, utr.stats.delta, **kwargs)
            header['rms'] = rms
            header['iter'] = np.nan
        else:
            raise ValueError('method must be \'iter\' or \'water\'')
        return cls(rf, header)
        
