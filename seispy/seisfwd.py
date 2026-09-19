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
            rftr = RFTrace.deconvolve(self.rstream[i], self.zstream[i], tshift=shift,
                                       f0=f0, **kwargs)
            rfstream.append(rftr)
        return rfstream

    def compute_vsapp(self, periods_s, *, method='iter', shift=10., f0=2.,
                      pre_filt=None, **kwargs):
        """Forward model apparent Vs for every stored P-wave ray parameter.

        :param periods_s: Increasing positive half-widths of the cosine-squared
            measurement window, in seconds.
        :param method: ``'iter'`` (default) or ``'water'`` deconvolution.
        :param shift: Time before the direct P arrival, in seconds.
        :param f0: Gaussian factor in the SeisPy deconvolution convention.
        :param pre_filt: Optional ``(freqmin, freqmax)`` bandpass in Hz, using
            two corners and zero phase, or ``None`` (default).
        :param kwargs: ``itmax``/``minderr`` for iterative deconvolution, or
            ``wlevel`` for water-level deconvolution.
        :return: Tuple of :class:`seispy.vsapp.Curve` in ``self.rayp`` order.

        Both R/Z and Z/Z are deconvolved with the same settings and amplitude
        scale. Waveforms are regenerated on a separate instance, so existing
        streams and the model are unchanged. Only ``ipha=1`` is supported.
        """
        from seispy.decon import deconit, deconwater
        from seispy.vsapp import _compute_vsapp_from_components

        self._check_vsapp_options(method, kwargs)
        if pre_filt is not None:
            pre_filt = np.asarray(pre_filt, dtype=float)
            if (pre_filt.shape != (2,) or not np.isfinite(pre_filt).all()
                    or not 0 < pre_filt[0] < pre_filt[1] < 0.5 / self.dt):
                raise ValueError('pre_filt must satisfy 0 < freqmin < freqmax < Nyquist')
        synthetic = SynSeis(self.depmod, self.rayp, self.dt, self.npts,
                            ipha=self.ipha)
        synthetic.run_fwd()
        if pre_filt is not None:
            synthetic.filter(*pre_filt, order=2, zerophase=True)
        curves = []
        for rayp, radial, vertical in zip(
                self.rayp, synthetic.rstream, synthetic.zstream, strict=True):
            if method == 'iter':
                rrf = deconit(radial.data, vertical.data, self.dt,
                              tshift=shift, f0=f0, phase='P', **kwargs)[0]
                zrf = deconit(vertical.data, vertical.data, self.dt,
                              tshift=shift, f0=f0, phase='P', **kwargs)[0]
            else:
                rrf = deconwater(radial.data, vertical.data, self.dt,
                                 tshift=shift, f0=f0, normalize=False,
                                 phase='P', **kwargs)[0]
                zrf = deconwater(vertical.data, vertical.data, self.dt,
                                 tshift=shift, f0=f0, normalize=False,
                                 phase='P', **kwargs)[0]
            times = np.arange(len(rrf)) * self.dt - shift
            curves.append(_compute_vsapp_from_components(times, rrf, zrf, rayp=rayp,
                                         periods_s=periods_s))
        return tuple(curves)

    def vsapp_kernel(self, periods_s, *, method='iter', shift=10., f0=2.,
                     pre_filt=None, zero_halfspace=False, **kwargs):
        """Return analytic apparent-Vs kernels for the stored ray parameters.

        Arguments match :meth:`compute_vsapp`; ``vp_vs_derivative`` and
        ``rho_vs_derivative`` may additionally supply local model slopes.
        ``periods_s`` contains half-widths of the cosine-squared measurement
        window in seconds.
        Results form a tuple in ``self.rayp`` order. Each Jacobian has one
        column per sampled model layer, including the half-space. Iterative
        gradients are local to the selected spikes and stopping iteration.
        Existing model arrays and synthetic streams are unchanged.

        :param zero_halfspace: Zero the last column of all returned Vs derivatives
            for inversion, retaining the full forward response, defaults to False
        :type zero_halfspace: bool, optional
        :return: One plottable VsappKernelResult per stored ray parameter
        :rtype: tuple[seispy.vsapp_kernel.VsappKernelResult]
        """
        self._check_vsapp_options(method, kwargs, kernel=True)
        return tuple(self.depmod.vsapp_kernel(
            rayp, periods_s, method=method, dt=self.dt, npts=self.npts,
            shift=shift, f0=f0, pre_filt=pre_filt, zero_halfspace=zero_halfspace, **kwargs,
        ) for rayp in self.rayp)

    def _check_vsapp_options(self, method, kwargs, kernel=False):
        if self.ipha != 1:
            raise ValueError('Vsapp supports only incident P waves (ipha=1)')
        if method not in ('iter', 'water'):
            raise ValueError("method must be 'iter' or 'water'")
        allowed = {'itmax', 'minderr'} if method == 'iter' else {'wlevel'}
        if kernel:
            allowed |= {'vp_vs_derivative', 'rho_vs_derivative'}
        unexpected = set(kwargs) - allowed
        if unexpected:
            raise TypeError('Unsupported Vsapp options: ' + ', '.join(sorted(unexpected)))
        try:
            arrays = [np.asarray(getattr(self.depmod, name), dtype=float)
                      for name in ('vp', 'vs', 'rho', 'thickness')]
        except (AttributeError, TypeError, ValueError) as exc:
            raise ValueError('Model must provide vp, vs, rho and thickness arrays') from exc
        if (any(a.ndim != 1 or not a.size or not np.isfinite(a).all() for a in arrays)
                or len({a.size for a in arrays}) != 1):
            raise ValueError('Model vp, vs, rho and thickness must be finite, nonempty, '
                             'one-dimensional arrays of equal length')
        vp, vs, rho, thickness = arrays
        if np.any(vs <= 0) or np.any(vp <= vs) or np.any(rho <= 0):
            raise ValueError('Model requires positive Vs and density, with Vp > Vs')
        if np.any(thickness[:-1] <= 0) or thickness[-1] != 0:
            raise ValueError('Finite layers require positive thickness and the last '
                             'half-space requires thickness=0')
        if not np.isscalar(self.dt) or not np.isfinite(self.dt) or self.dt <= 0:
            raise ValueError('dt must be finite and positive')
        if (isinstance(self.npts, (bool, np.bool_))
                or not isinstance(self.npts, (int, np.integer)) or self.npts < 8):
            raise ValueError('npts must be an integer >= 8')
        rayp = np.asarray(self.rayp, dtype=float)
        if (rayp.ndim != 1 or not rayp.size or not np.isfinite(rayp).all()
                or np.any(rayp <= 0)):
            raise ValueError('rayp must be a nonempty one-dimensional array of '
                             'finite positive slownesses in s/km')
        if np.any(rayp[:, None] * vp[None, :] >= 1):
            raise ValueError('Vsapp supports only subcritical incidence (rayp*Vp < 1)')
