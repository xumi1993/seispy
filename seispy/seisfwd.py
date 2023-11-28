import numpy as np
from scipy.fftpack import ifft
from obspy.signal.util import next_pow_2
from seispy.utils import scalar_instance, array_instance
from obspy import Trace, Stream
from seispy.decon import RFTrace
from numba import njit

ei = 0+1j


@njit(fastmath=True,cache=True)
def e_inverse(omega, rho, alpha, beta, p):
    """ E_inverse (Aki & Richards, pp. 161, Eq. (5.71))

    Parameters
    ----------
    omega : _type_
        _description_
    rho : _type_
        _description_
    alpha : _type_
        _description_
    beta : _type_
        _description_
    p : _type_
        _description_
    """
    e_inv = np.zeros((4,4), dtype=np.complex128)
    eta = np.sqrt(1.0/(beta*beta) - p*p)
    xi  = np.sqrt(1.0/(alpha*alpha) - p*p)
    bp = 1.0 - 2.0*beta*beta*p*p

    e_inv[0,0] = beta*beta*p/alpha
    e_inv[0,1] = bp/(2.0*alpha*xi)
    e_inv[0,2] = -p/(2.0*omega*rho*alpha*xi) * ei
    e_inv[0,3] = -1.0/(2.0*omega*rho*alpha) * ei
    e_inv[1,0] = bp / (2.0*beta*eta)
    e_inv[1,1] = -beta*p
    e_inv[1,2] = -1.0/(2.0*omega*rho*beta) * ei
    e_inv[1,3] = p/(2.0*omega*rho*beta*eta) * ei
    e_inv[2,0] = e_inv[0,0]
    e_inv[2,1] = - e_inv[0,1]
    e_inv[2,2] = - e_inv[0,2]
    e_inv[2,3] = e_inv[0,3]
    e_inv[3,0] = e_inv[1,0]
    e_inv[3,1] = - e_inv[1,1]
    e_inv[3,2] = - e_inv[1,2]
    e_inv[3,3] = e_inv[1,3]
    return e_inv

@njit(fastmath=True,cache=True)
def propagator_sol(omega, rho, alpha, beta, p, z):
    """ 
    propagator (Aki & Richards, pp. 398, Eq. (3) in Box 9.1)

    Parameters
    ----------
    omega : _type_
        _description_
    rho : _type_
        _description_
    alpha : _type_
        _description_
    beta : _type_
        _description_
    p : _type_
        _description_
    z : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    
    p_mat = np.zeros((4,4), dtype=np.complex128)
    beta2 = beta*beta
    p2 = p*p
    bp = 1.0 -2.0*beta2*p2
    eta = np.sqrt(1.0/(beta2) - p2)
    xi  = np.sqrt(1.0/(alpha*alpha) - p2)
    cos_xi = np.cos(omega*xi*z)
    cos_eta = np.cos(omega*eta*z)
    sin_xi = np.sin(omega*xi*z)
    sin_eta = np.sin(omega*eta*z)

    p_mat[0,0] = 2.0*beta2*p2*cos_xi + bp*cos_eta
    p_mat[0,1] = p*( bp/xi*sin_xi - 2.0*beta2*eta*sin_eta ) * ei
    p_mat[0,2] = (p2/xi*sin_xi + eta*sin_eta)/(omega*rho)
    p_mat[0,3] = p*(-cos_xi + cos_eta)/(omega*rho) * ei  
    p_mat[1,0] = p*( 2.0*beta2*xi*sin_xi - bp/eta*sin_eta ) * ei
    p_mat[1,1] = bp*cos_xi + 2.0*beta2*p2*cos_eta
    p_mat[1,2] = p_mat[0,3]
    p_mat[1,3] = (xi*sin_xi + p2/eta*sin_eta)/(omega*rho)
    p_mat[2,0] = omega*rho*( -4.0*beta2*beta2*p2*xi*sin_xi - bp*bp/eta*sin_eta )
    p_mat[2,1] = 2.0*omega*beta2*rho*p*bp*( cos_xi - cos_eta ) * ei
    p_mat[2,2] = p_mat[0,0]
    p_mat[2,3] = p_mat[1,0]
    p_mat[3,0] = p_mat[2,1]
    p_mat[3,1] = -omega*rho*( bp*bp/xi*sin_xi + 4.0*beta2*beta2*p2*eta*sin_eta  )
    p_mat[3,2] = p_mat[0,1]  
    p_mat[3,3] = p_mat[1,1]

    return p_mat

@njit(fastmath=True,cache=True)
def haskell(omega, p, nl, ipha, alpha, beta, rho, h):
    i0 = 0
    e_inv = e_inverse(omega, rho[-1], alpha[-1], beta[-1], p)
    p_mat = propagator_sol(omega, rho[i0], alpha[i0], beta[i0], p, h[i0] )
    for i in range(i0+1, nl):
        p_mat2 = propagator_sol(omega, rho[i], alpha[i], beta[i], p, h[i])
        p_mat = p_mat2 @ p_mat
    if nl > i0+1:
        sl = e_inv @ p_mat
    else:
        sl = e_inv
    denom = sl[2,0] * sl[3,1] - sl[2,1] * sl[3,0]
    if ipha >= 0:
        ur = sl[3,1] / denom
        uz = - sl[3,0] / denom
    else:
        ur = - sl[2,1] / denom
        uz = sl[2,0] / denom
    return ur, uz

@njit(fastmath=True,cache=True)
def fwd_seis(rayp, dt, npts, ipha, alpha, beta, rho, h):
    nlay = h.size
    ur_freq = np.zeros(npts, dtype=np.complex128)
    uz_freq = np.zeros(npts, dtype=np.complex128)
    nhalf = int(npts / 2 + 1)
    for i in range(1, nhalf):
        omg = 2*np.pi * i / (npts * dt)
        ur_freq[i], uz_freq[i] = haskell(omg, rayp, nlay, ipha, 
                                         alpha, beta, rho, h)
    return ur_freq, uz_freq


class SynSeis():
    def __init__(self, depmod, rayp, dt, npts=2500, ipha=1, filter=None) -> None:
        """_summary_

        Parameters
        ----------
        depmod : _type_
            DepModel class
        rayp : _type_
            Ray-parameter in s/km
        dt : _type_
            Time interval
        npts : _type_
            samples of synthetic waveform
        ipha : _type_
            Specify incident wave 1 for P and -1 for S
        """
        self.depmod = depmod
        self.dt = dt
        self.npts = npts
        if not (array_instance(rayp) or scalar_instance(rayp)):
            raise TypeError('The rayp should be in float, list and np.ndarray')
        if scalar_instance(rayp):
            self.rayp = [rayp]
        else:
            self.rayp = rayp
        self.ipha = ipha

    def run_fwd(self):
        """Forward modelling synthetic seismograms.

        ``SynSeis.rstream`` and ``SynSeis.zstream`` are generated as 
        radial and vertical Seismograms in ``Obspy.Stream`` type.
        """
        self.rstream = Stream()
        self.zstream = Stream()
        npts_max = next_pow_2(self.npts)
        for _, rayp in enumerate(self.rayp):
            ur_freq, uz_freq = fwd_seis(rayp, self.dt, npts_max, self.ipha,
                            self.depmod.vp, self.depmod.vs, self.depmod.rho,
                            self.depmod.thickness)
            ur = ifft(ur_freq).real[::-1]/npts_max
            uz = -ifft(uz_freq).real[::-1]/npts_max
            tr = Trace(data=ur)
            tr.stats.delta = self.dt
            self.rstream.append(tr)
            tr = Trace(data=uz)
            tr.stats.delta = self.dt
            self.zstream.append(tr)

    def filter(self, freqmin, freqmax, order=2, zerophase=True):
        """Apply a bandpass filter on synthetic waveforms

        Parameters
        ----------
        freqmin : float
            Minimum cut-off frequency
        freqmax : float
            maximum cut-off frequency
        order : int, optional
            Order of filter, by default 2
        zerophase : bool, optional
            whether use a zero-phase filter, by default True
        """
        for st in [self.rstream, self.zstream]:
            st.filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                      corners=order, zerophase=zerophase)

    def run_deconvolution(self, pre_filt=[0.05, 2], shift=10, f0=2.0, **kwargs):
        if pre_filt is not None:
            self.filter(*pre_filt)
        rfstream = Stream()
        for i, _ in enumerate(self.rayp):
            rftr = RFTrace.deconvolute(self.rstream[i], self.zstream[i], tshift=shift,
                                       f0=f0, **kwargs)
            rfstream.append(rftr)
        return rfstream

