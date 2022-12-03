import numpy as np


def smooth(x, half_len=5, window='flat'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        helf_len: the half dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    window_len = 2*half_len+1

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming',"
                         "'bartlett', 'blackman'")
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':
        w = np.ones(window_len, 'd')  # moving average
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[half_len:-half_len]


def whiten(data, Nfft, delta, f1, f2, f3, f4):
    """This function takes 1-dimensional *data* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *freqmin* and *freqmax*
    and returns the whitened fft.

    :type data: :class:`numpy.ndarray`
    :param data: Contains the 1D time series to whiten
    :type Nfft: int
    :param Nfft: The number of points to compute the FFT
    :type delta: float
    :param delta: The sampling frequency of the `data`
    :type freqmin: float
    :param freqmin: The lower frequency bound
    :type freqmax: float
    :param freqmax: The upper frequency bound
    :type plot: bool
    :param plot: Whether to show a raw plot of the action (default: False)

    :rtype: :class:`numpy.ndarray`
    :returns: The FFT of the input trace, whitened between the frequency bounds
"""
    dom = 1/delta/Nfft
    nt1 = int(f1/dom)
    nt2 = int(f2/dom)
    nt3 = int(f3/dom)
    nt4 = int(f4/dom)

    FFTRawSign = np.fft.fft(data, Nfft)

    FFTRawSign /= smooth(np.abs(FFTRawSign), half_len=20)
    # Left tapering:
    FFTRawSign[0:nt1] *= 0
    FFTRawSign[nt1:nt2] = np.cos(np.linspace(np.pi / 2., np.pi, nt2 - nt1)) ** 2 * np.exp(1j * np.angle(FFTRawSign[nt1:nt2]))
    #FFTRawSign[nt1:nt2] = np.cos(np.linspace(np.pi / 2., np.pi, nt2 - nt1)) ** 2 * FFTRawSign[nt1:nt2]
    # Pass band:
    FFTRawSign[nt2:nt3] = np.exp(1j * np.angle(FFTRawSign[nt2:nt3]))
    # Right tapering:
    FFTRawSign[nt3:nt4] = np.cos(np.linspace(0., np.pi / 2., nt4 - nt3)) ** 2 * np.exp(1j * np.angle(FFTRawSign[nt3:nt4]))
    #FFTRawSign[nt3:nt4] = np.cos(np.linspace(0., np.pi / 2., nt4 - nt3)) ** 2 * FFTRawSign[nt3:nt4]
    FFTRawSign[nt4:Nfft+1] *= 0
    
    # Hermitian symmetry (because the input is real)
    FFTRawSign[int(-Nfft/2+1):] = FFTRawSign[1:int(Nfft/2)].conjugate()[::-1]
    
    return FFTRawSign, dom

