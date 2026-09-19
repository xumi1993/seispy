"""Analytic layered-model RF and apparent-Vs Jacobians.

Use vsapp_kernel(DepModel, rayp, periods_s) for SeisPy models, or the
array-level forward_and_jacobian_iter / forward_and_jacobian functions.
Matrix derivatives use cached prefix/suffix products, and iterative
filtering/correlation/phase conventions reuse seispy.decon directly.

Derivatives are local to the active water-level or iterative execution
branch. Iterative peak selection and stopping may change under finite
perturbations. branch_warning=False does not certify branch stability.

Layer columns refer to the current DepModel.vs grid, including its final
zero-thickness half-space. By default Vp, density and thickness are fixed;
explicit local Vp/rho slopes allow constrained total derivatives. No RFapp
installation, finite differences, or external benchmark code is required.
"""

from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.fftpack import fft, ifft

from seispy.decon import correl as _correl
from seispy.decon import gaussFilter as _gauss_filter
from seispy.decon import gfilter as _gfilter
from seispy.decon import phaseshift as _phaseshift
from seispy.vsapp import _compute_vsapp_from_components

__all__ = [
    "vsapp_kernel",
    "VsappKernelResult",
    "IterDeconResult",
    "forward_and_jacobian_iter",
    "forward_and_jacobian",
    "deconit_tangent",
]


@dataclass(frozen=True)
class VsappKernelResult:
    """Synthetic RFs, apparent Vs and its layer sensitivities.

    ``jacobian[i, j]`` is dVsapp(T_i)/dVs_j, with shape (periods, layers)
    and units (km/s)/(km/s). Columns include the final half-space. These
    are derivatives with respect to whole-layer velocities, not sensitivity
    per unit depth. ``periods_s`` and ``thickness_km`` retain independent
    copies of the forward model's window half-widths and layer thicknesses.
    The last thickness is zero and denotes the infinite half-space.
    With ``diagnostics['zero_halfspace']`` enabled, the last columns of
    ``jacobian``, ``dradial_dvs`` and ``dvertical_dvs`` are zeroed for inversion.
    The synthetic RFs and apparent velocities retain the full model response.
    """

    times: np.ndarray
    radial: np.ndarray
    vertical: np.ndarray
    vsapp: np.ndarray
    dradial_dvs: np.ndarray
    dvertical_dvs: np.ndarray
    jacobian: np.ndarray
    diagnostics: dict
    periods_s: np.ndarray
    thickness_km: np.ndarray

    def plot(self, ax=None, *, cmap='RdBu_r', vmax=None, max_depth=None,
             colorbar=True, title=None, show=False):
        """Plot the Vsapp sensitivity colormap against period and depth.

        Depth increases downward from the top of the forward model. Each
        cell shows the unscaled derivative for a whole model layer; values
        are not divided by layer thickness or normalized by period. Period
        cell edges are midpoints between the supplied window half-widths.
        The color scale is symmetric about zero.

        The final half-space is shown down to max_depth and labelled as
        such. Its displayed extent is only a plotting convention, not an
        additional finite model layer. By default it occupies the thickness
        of the preceding layer, or 1 km for a homogeneous half-space.

        :param ax: Axes to draw on. Create a new figure if None, defaults to None
        :type ax: matplotlib.axes.Axes, optional
        :param cmap: Matplotlib colormap, defaults to 'RdBu_r'
        :type cmap: str or matplotlib.colors.Colormap, optional
        :param vmax: Positive absolute color limit. None uses the largest
            absolute sensitivity; colors span [-vmax, vmax], defaults to None
        :type vmax: float, optional
        :param max_depth: Maximum displayed depth in km below the model top.
            None includes all finite layers and a half-space band, defaults to None
        :type max_depth: float, optional
        :param colorbar: Add a labelled colorbar, defaults to True
        :type colorbar: bool, optional
        :param title: Plot title, defaults to None
        :type title: str, optional
        :param show: Show the figure after drawing, defaults to False
        :type show: bool, optional
        :return: Figure and axes for further styling or saving with fig.savefig
        :rtype: (matplotlib.figure.Figure, matplotlib.axes.Axes)
        :raises ValueError: If coordinates, sensitivities or plot limits are invalid
        """
        import matplotlib.pyplot as plt

        periods = _array(self.periods_s, 'periods_s')
        thickness = _array(self.thickness_km, 'thickness_km')
        if np.any(periods <= 0) or np.any(np.diff(periods) <= 0):
            raise ValueError('periods_s must be positive and strictly increasing')
        if np.any(thickness[:-1] <= 0) or thickness[-1] != 0:
            raise ValueError('Finite layers need positive thickness; '
                             'half-space thickness must be zero')
        values = np.asarray(self.jacobian, dtype=float)
        if values.shape != (periods.size, thickness.size) or not np.isfinite(values).all():
            raise ValueError('jacobian must be finite with shape (periods, layers)')
        limit = float(np.max(np.abs(values))) or 1.0
        if vmax is not None:
            limit = float(vmax)
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError('vmax must be finite and positive')

        layer_tops = np.r_[0., np.cumsum(thickness[:-1])]
        halfspace_top = layer_tops[-1]
        default_bottom = halfspace_top + (thickness[-2] if thickness.size > 1 else 1.)
        bottom = default_bottom if max_depth is None else float(max_depth)
        if not np.isfinite(bottom) or bottom <= 0:
            raise ValueError('max_depth must be finite and positive')
        depth_edges = np.r_[layer_tops, max(bottom, default_bottom)]
        if periods.size == 1:
            period_edges = np.array([periods[0] / 2, periods[0] * 1.5])
        else:
            period_edges = np.r_[max(0., periods[0] - (periods[1] - periods[0]) / 2),
                                 (periods[:-1] + periods[1:]) / 2,
                                 periods[-1] + (periods[-1] - periods[-2]) / 2]

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
        else:
            fig = ax.figure
        mesh = ax.pcolormesh(period_edges, depth_edges, values.T, shading='flat',
                             cmap=cmap, vmin=-limit, vmax=limit, rasterized=True)
        ax.set_xlabel('Period (s)')
        ax.set_ylabel('Depth below model top (km)')
        ax.set_xlim(period_edges[0], period_edges[-1])
        ax.set_ylim(bottom, 0.)
        if bottom > halfspace_top:
            if halfspace_top > 0:
                ax.axhline(halfspace_top, color='0.35', linewidth=0.7, linestyle='--')
            ax.text(.98, (halfspace_top + bottom) / 2, 'Half-space',
                    transform=ax.get_yaxis_transform(), ha='right', va='center',
                    fontsize='small', bbox=dict(facecolor='white', alpha=.8, edgecolor='none'))
        if colorbar:
            fig.colorbar(mesh, ax=ax, label=r'$\partial V_{S,\mathrm{app}}/\partial V_{S,j}$')
        if title is not None:
            ax.set_title(title)
        if show:
            plt.show()
        return fig, ax


@njit(cache=True)
def _propagator_tangent(omega, rho, alpha, beta, p, h, da, dr):
    """SeisPy propagator and dP/dVs with dBeta=1, dAlpha=da, dRho=dr."""
    u = beta * beta
    du = 2.0 * beta
    p2 = p * p
    b = 1.0 - 2.0 * u * p2
    db = -2.0 * du * p2
    eta = np.sqrt(1.0 / u - p2)
    xi = np.sqrt(1.0 / (alpha * alpha) - p2)
    de = -1.0 / (beta**3 * eta)
    dx = -da / (alpha**3 * xi)
    cx, ce = np.cos(omega * xi * h), np.cos(omega * eta * h)
    sx, se = np.sin(omega * xi * h), np.sin(omega * eta * h)
    dcx, dce = -sx * omega * h * dx, -se * omega * h * de
    dsx, dse = cx * omega * h * dx, ce * omega * h * de
    x1, x2, e1, e2 = sx / xi, xi * sx, se / eta, eta * se
    dx1, dx2 = dsx / xi - sx * dx / xi**2, dx * sx + xi * dsx
    de1, de2 = dse / eta - se * de / eta**2, de * se + eta * dse
    q, dq = omega * rho, omega * dr
    c, dc = cx - ce, dcx - dce
    a = np.zeros((4, 4), dtype=np.complex128)
    d = np.zeros((4, 4), dtype=np.complex128)

    a[0, 0] = 2 * u * p2 * cx + b * ce
    d[0, 0] = 2 * p2 * (du * cx + u * dcx) + db * ce + b * dce
    a[0, 1] = 1j * p * (b * x1 - 2 * u * e2)
    d[0, 1] = 1j * p * (db * x1 + b * dx1 - 2 * (du * e2 + u * de2))
    a[0, 2] = (p2 * x1 + e2) / q
    d[0, 2] = (p2 * dx1 + de2) / q - a[0, 2] * dq / q
    a[0, 3] = -1j * p * c / q
    d[0, 3] = -1j * p * dc / q - a[0, 3] * dq / q
    a[1, 0] = 1j * p * (2 * u * x2 - b * e1)
    d[1, 0] = 1j * p * (2 * (du * x2 + u * dx2) - db * e1 - b * de1)
    a[1, 1] = b * cx + 2 * u * p2 * ce
    d[1, 1] = db * cx + b * dcx + 2 * p2 * (du * ce + u * dce)
    a[1, 2], d[1, 2] = a[0, 3], d[0, 3]
    a[1, 3] = (x2 + p2 * e1) / q
    d[1, 3] = (dx2 + p2 * de1) / q - a[1, 3] * dq / q
    a[2, 0] = -q * (4 * u * u * p2 * x2 + b * b * e1)
    d[2, 0] = a[2, 0] * dq / q - q * (
        4 * p2 * (2 * u * du * x2 + u * u * dx2) + 2 * b * db * e1 + b * b * de1
    )
    a[2, 1] = 2j * q * u * p * b * c
    d[2, 1] = 2j * p * (
        dq * u * b * c + q * du * b * c + q * u * db * c + q * u * b * dc
    )
    a[2, 2], d[2, 2] = a[0, 0], d[0, 0]
    a[2, 3], d[2, 3] = a[1, 0], d[1, 0]
    a[3, 0], d[3, 0] = a[2, 1], d[2, 1]
    a[3, 1] = -q * (b * b * x1 + 4 * u * u * p2 * e2)
    d[3, 1] = a[3, 1] * dq / q - q * (
        2 * b * db * x1 + b * b * dx1 + 4 * p2 * (2 * u * du * e2 + u * u * de2)
    )
    a[3, 2], d[3, 2] = a[0, 1], d[0, 1]
    a[3, 3], d[3, 3] = a[1, 1], d[1, 1]
    return a, d


@njit(cache=True)
def _e_inverse_tangent(omega, rho, alpha, beta, p, da, dr):
    """Half-space boundary matrix and its local constrained Vs derivative."""
    eta = np.sqrt(1.0 / beta**2 - p * p)
    xi = np.sqrt(1.0 / alpha**2 - p * p)
    de, dx = -1.0 / (beta**3 * eta), -da / (alpha**3 * xi)
    b, db = 1.0 - 2 * beta**2 * p * p, -4 * beta * p * p
    e = np.zeros((4, 4), dtype=np.complex128)
    d = np.zeros((4, 4), dtype=np.complex128)
    e[0, 0] = beta**2 * p / alpha
    d[0, 0] = 2 * beta * p / alpha - e[0, 0] * da / alpha
    e[0, 1] = b / (2 * alpha * xi)
    d[0, 1] = db / (2 * alpha * xi) - e[0, 1] * (da / alpha + dx / xi)
    e[0, 2] = -1j * p / (2 * omega * rho * alpha * xi)
    d[0, 2] = -e[0, 2] * (dr / rho + da / alpha + dx / xi)
    e[0, 3] = -1j / (2 * omega * rho * alpha)
    d[0, 3] = -e[0, 3] * (dr / rho + da / alpha)
    e[1, 0] = b / (2 * beta * eta)
    d[1, 0] = db / (2 * beta * eta) - e[1, 0] * (1 / beta + de / eta)
    e[1, 1], d[1, 1] = -beta * p, -p
    e[1, 2] = -1j / (2 * omega * rho * beta)
    d[1, 2] = -e[1, 2] * (dr / rho + 1 / beta)
    e[1, 3] = 1j * p / (2 * omega * rho * beta * eta)
    d[1, 3] = -e[1, 3] * (dr / rho + 1 / beta + de / eta)
    for j in range(4):
        sign = -1.0 if j == 1 or j == 2 else 1.0
        e[2, j], d[2, j] = sign * e[0, j], sign * d[0, j]
        e[3, j], d[3, j] = sign * e[1, j], sign * d[1, j]
    return e, d


@njit(cache=True)
def _spectra_and_tangents(vp, vs, rho, h, da, dr, p, dt, nfft):
    """O(nfft*nlayer) complete forward and all Vs columns via cached products."""
    nl = len(vs)
    ur, uz = np.zeros(nfft, np.complex128), np.zeros(nfft, np.complex128)
    dur = np.zeros((nfft, nl), np.complex128)
    duz = np.zeros((nfft, nl), np.complex128)
    for k in range(1, nfft // 2 + 1):
        omega = 2 * np.pi * k / (nfft * dt)
        e, de = _e_inverse_tangent(omega, rho[-1], vp[-1], vs[-1], p, da[-1], dr[-1])
        matrices = np.empty((nl, 4, 4), np.complex128)
        derivs = np.empty((nl, 4, 4), np.complex128)
        prefix = np.empty((nl + 1, 4, 4), np.complex128)
        prefix[0] = np.eye(4, dtype=np.complex128)
        for j in range(nl):
            matrices[j], derivs[j] = _propagator_tangent(
                omega, rho[j], vp[j], vs[j], p, h[j], da[j], dr[j]
            )
            prefix[j + 1] = matrices[j] @ prefix[j]
        s = e @ prefix[nl]
        det = s[2, 0] * s[3, 1] - s[2, 1] * s[3, 0]
        ur[k], uz[k] = s[3, 1] / det, -s[3, 0] / det
        left = e
        for j in range(nl - 1, -1, -1):
            ds = left @ derivs[j] @ prefix[j]
            if j == nl - 1:
                ds += de @ prefix[nl]
            ddet = (
                ds[2, 0] * s[3, 1] + s[2, 0] * ds[3, 1]
                - ds[2, 1] * s[3, 0] - s[2, 1] * ds[3, 0]
            )
            dur[k, j] = (ds[3, 1] - ur[k] * ddet) / det
            duz[k, j] = (-ds[3, 0] - uz[k] * ddet) / det
            left = left @ matrices[j]
    return ur, uz, dur, duz


def _water_level_tangents(r, z, dr, dz, dt, shift, f0, wlevel):
    """Differentiate both R/Z and Z/Z, including the moving water-level floor."""
    n = len(r)
    u, w = fft(r), fft(z)
    du, dw = fft(dr, axis=0), fft(dz, axis=0)
    power = (w * w.conjugate()).real
    dpower = 2 * (w.conjugate()[:, None] * dw).real
    peak_index = int(np.argmax(power))
    peak = power[peak_index]
    if peak <= 0 or not np.isfinite(peak):
        raise ValueError("Vertical waveform spectrum is zero or nonfinite")
    floor = wlevel * peak
    active = power < floor
    denominator = np.maximum(power, floor)
    ddenominator = np.where(active[:, None], wlevel * dpower[peak_index], dpower)
    omega = 2 * np.pi * np.fft.fftfreq(n, dt)
    omega[n // 2] = abs(omega[n // 2])  # SeisPy uses positive Nyquist for its shift.
    multiplier = np.exp(-0.25 * (omega / f0) ** 2) / dt * np.exp(-1j * omega * shift)

    numerator = u * w.conjugate()
    dnumerator = du * w.conjugate()[:, None] + u[:, None] * dw.conjugate()
    ratio = numerator / denominator
    dratio = (dnumerator - ratio[:, None] * ddenominator) / denominator[:, None]
    zratio = power / denominator
    dzratio = (dpower - zratio[:, None] * ddenominator) / denominator[:, None]
    radial = ifft(multiplier * ratio).real
    vertical = ifft(multiplier * zratio).real
    dradial = ifft(multiplier[:, None] * dratio, axis=0).real
    dvertical = ifft(multiplier[:, None] * dzratio, axis=0).real

    # +/- frequency ties are expected for real inputs and have identical slopes.
    ties = np.flatnonzero(np.abs(power - peak) <= 16 * np.finfo(float).eps * peak)
    tie_slope_spread = float(np.max(np.abs(dpower[ties] - dpower[peak_index])) / peak)
    switch_margin = float(np.min(np.abs(power - floor)) / peak)
    diagnostics = {
        "waterlevel_active_bins": int(active.sum()),
        "waterlevel_switch_margin_relative": switch_margin,
        "spectral_max_tie_count": int(len(ties)),
        "spectral_max_tie_slope_spread_relative": tie_slope_spread,
        "branch_warning": bool(switch_margin < 1e-12 or tie_slope_spread > 1e-8),
    }
    return radial, vertical, dradial, dvertical, diagnostics


def _integration_operator(times, periods):
    """Exact linear weights of the paired-component interpolated trapezoids."""
    operator = np.zeros((len(periods), len(times)))
    for i, period in enumerate(periods):
        inside = times[(times > -period) & (times < period)]
        grid = np.unique(np.concatenate(([-period, 0.0, period], inside)))
        widths = np.diff(grid)
        quadrature = np.zeros(len(grid))
        quadrature[:-1] += widths / 2
        quadrature[1:] += widths / 2
        quadrature *= np.cos(np.pi * grid / (2 * period)) ** 2
        quadrature[[0, -1]] = 0
        right = np.clip(np.searchsorted(times, grid), 1, len(times) - 1)
        fraction = (grid - times[right - 1]) / (times[right] - times[right - 1])
        np.add.at(operator[i], right - 1, quadrature * (1 - fraction))
        np.add.at(operator[i], right, quadrature * fraction)
    return operator


def _array(value, name, length=None):
    array = np.array(value, dtype=float, copy=True)
    if array.ndim != 1 or not len(array) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a nonempty finite 1-D array")
    if length is not None and len(array) != length:
        raise ValueError(f"{name} must have {length} entries")
    return array


def forward_and_jacobian(
    vp, vs, rho, thickness, *, p=0.06, dt=0.05, npts=1024, shift=10.0,
    f0=2.0, wlevel=0.05, periods=None, vp_vs_derivative=None, rho_vs_derivative=None,
    pre_filt=None, zero_halfspace=False,
):
    """Return a water-level synthetic RF, Vsapp, and exact local Vs Jacobians.

    Units: km/s, g/cm^3, km, s/km, s. Slopes default to zero (fixed Vp/rho).
    pre_filt optionally applies SynSeis's second-order zero-phase bandpass
    to waveforms and every tangent before deconvolution.
    A supplied slope array is dVp_j/dVs_j or dRho_j/dVs_j. Cross-layer
    constraints and layer-thickness derivatives are outside this module.
    Times are relative to P; the actual length follows SeisPy's next_pow_2.
    The derivative is valid on the selected water-level/max branches; inspect
    diagnostics['branch_warning'] before interpreting it at a branch boundary.

    :param zero_halfspace: Zero the half-space column of all returned Vs
        derivatives without changing the forward response, defaults to False
    :type zero_halfspace: bool, optional
    """
    if not isinstance(zero_halfspace, (bool, np.bool_)):
        raise TypeError('zero_halfspace must be a bool')
    vs = _array(vs, "vs")
    nl = len(vs)
    vp, rho, h = (_array(a, name, nl) for a, name in
                  [(vp, "vp"), (rho, "rho"), (thickness, "thickness")])
    if np.any(vs <= 0) or np.any(vp <= vs) or np.any(rho <= 0):
        raise ValueError("Require positive Vs/rho and Vp > Vs")
    if np.any(h[:-1] <= 0) or h[-1] != 0:
        raise ValueError("Finite layers need positive thickness; last half-space needs thickness=0")
    for value, name in [(p, "p"), (dt, "dt"), (f0, "f0"), (wlevel, "wlevel")]:
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if wlevel >= 1:
        raise ValueError("wlevel must be < 1")
    if np.any(vp * p >= 1) or np.any(vs * p >= 1):
        raise ValueError("Only subcritical incidence p*Vp < 1 and p*Vs < 1 is supported")
    if isinstance(npts, (bool, np.bool_)) or not isinstance(npts, (int, np.integer)) or npts < 8:
        raise ValueError("npts must be an integer >= 8")
    if not np.isfinite(shift) or shift <= 0:
        raise ValueError("shift must be finite and positive")
    if f0 * dt > 0.5:
        raise ValueError("Gaussian pulse is undersampled (f0*dt > 0.5)")
    periods = _array(np.linspace(0.1, 3.5, 35) if periods is None else periods, "periods")
    if np.any(periods <= 0) or np.any(np.diff(periods) <= 0) or periods[0] < 2 * dt - 1e-9:
        raise ValueError("periods must increase strictly and be >= 2*dt")
    nfft = 1 << (int(npts) - 1).bit_length()
    times = np.arange(nfft) * dt - shift
    if times[0] > -periods[-1] + 1e-8 or times[-1] < periods[-1] - 1e-8:
        raise ValueError("Waveform time span must cover every [-T,T] window")
    da = np.zeros(nl) if vp_vs_derivative is None else _array(vp_vs_derivative, "dVp/dVs", nl)
    dr = np.zeros(nl) if rho_vs_derivative is None else _array(rho_vs_derivative, "dRho/dVs", nl)
    ur, uz, dur, duz = _spectra_and_tangents(vp, vs, rho, h, da, dr, p, dt, nfft)
    # Preserve SeisPy's one-sided spectra, sign, real projection, time reversal,
    # and extra normalization; blindly using irfft would define another forward.
    r, z = ifft(ur).real[::-1] / nfft, -ifft(uz).real[::-1] / nfft
    d_r, d_z = ifft(dur, axis=0).real[::-1] / nfft, -ifft(duz, axis=0).real[::-1] / nfft
    r, z, d_r, d_z = _prefilter_waveforms(r, z, d_r, d_z, dt, pre_filt)
    radial, vertical, d_radial, d_vertical, diagnostics = _water_level_tangents(
        r, z, d_r, d_z, dt, shift, f0, wlevel
    )
    curve = _compute_vsapp_from_components(times, radial, vertical, rayp=p, periods_s=periods)
    if any(status != "ok" for status in curve.status):
        raise ValueError(f"Vsapp is outside the valid measurement branch: {curve.status}")
    operator = _integration_operator(times, periods)
    a, b = curve.radial_integral, curve.vertical_integral
    # A direct contraction also handles an identically zero Z/Z tangent. Avoid
    # the platform BLAS matmul path, which emits spurious floating exceptions
    # for the (nfft, 1) half-space case on the tested Apple/NumPy runtime.
    da_int = np.einsum("it,tj->ij", operator, d_radial, optimize=False)
    db_int = np.einsum("it,tj->ij", operator, d_vertical, optimize=False)
    prefactor = np.cos(0.5 * np.arctan(a / b)) / (2 * p * (a * a + b * b))
    jacobian = prefactor[:, None] * (b[:, None] * da_int - a[:, None] * db_int)
    if not all(np.isfinite(x).all() for x in (radial, vertical, d_radial, d_vertical, jacobian)):
        raise ValueError("Nonfinite forward or tangent result; model is numerically singular")
    diagnostics.update({"nfft": nfft, "n_layers_including_halfspace": nl,
                        "parameterization": "Vs with specified local Vp/rho slopes",
                        "derivative_method": "analytic matrix derivatives and product rule",
                        "deconvolution": "water",
                        "prefilter": None if pre_filt is None else list(pre_filt)})
    diagnostics['zero_halfspace'] = bool(zero_halfspace)
    if zero_halfspace:
        jacobian[:, -1] = 0.
        d_radial[:, -1] = 0.
        d_vertical[:, -1] = 0.
    return VsappKernelResult(times, radial, vertical, curve.vs_km_s, d_radial, d_vertical,
                             jacobian, diagnostics, periods.copy(), h.copy())


@dataclass(frozen=True)
class IterDeconResult:
    rf: np.ndarray
    jacobian: np.ndarray
    path: list[int]
    rms: np.ndarray
    improvements: np.ndarray
    peak_margins: np.ndarray
    amplitudes: np.ndarray
    stop_margins: np.ndarray
    degenerate_tangent_iterations: list[int]
    forced: bool


def _columns_filter(x, spectrum, dt):
    return ifft(fft(x, axis=0) * spectrum[:, None] * dt, axis=0).real


def deconit_tangent(
    uin, win, duin, dwin, dt, *, tshift=10.0, f0=2.0, itmax=400,
    minderr=0.001, forced_path=None,
):
    """Analytic fixed-branch tangent of deconit(..., phase='P', nt=None).

    uin/win shape (nt,); duin/dwin shape (nt, n_parameters). The returned rms
    retains EVERY executed iteration, unlike SeisPy's legacy truncated rms.
    path length is the iteration count (SeisPy returns last zero-based index).
    amplitudes are normalized correlations c_l = dt * spike_amplitude.
    forced_path is a diagnostic only: replay exactly those lags and that count,
    bypassing argmax and stopping. The normal mode uses the real SeisPy rules.
    """
    u = _array(uin, "uin")
    w = _array(win, "win", len(u))
    nt = len(u)
    nfft = 1 << (nt - 1).bit_length()
    if nt < 8:
        raise ValueError("At least 8 waveform samples are required")
    for value, name in [(dt, "dt"), (f0, "f0")]:
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not np.isfinite(tshift) or tshift < 0:
        raise ValueError("tshift must be finite and nonnegative")
    if not np.isfinite(minderr) or minderr < 0:
        raise ValueError("minderr must be finite and nonnegative")
    if isinstance(itmax, (bool, np.bool_)) or not isinstance(itmax, (int, np.integer)) or itmax < 1:
        raise ValueError("itmax must be a positive integer")
    du = np.array(duin, dtype=float, copy=True)
    dw = np.array(dwin, dtype=float, copy=True)
    if du.ndim != 2 or du.shape[0] != nt or dw.shape != du.shape:
        raise ValueError("duin and dwin must have identical (nt, n_parameters) shapes")
    if not np.isfinite(du).all() or not np.isfinite(dw).all():
        raise ValueError("Waveform derivatives must be finite")
    npart = du.shape[1]
    maxlag = nfft // 2 - 1  # exact P search slice in deconit
    if forced_path is not None:
        forced_path = list(forced_path)
        if not forced_path or any(
            isinstance(k, (bool, np.bool_)) or not isinstance(k, (int, np.integer))
            or k < 0 or k >= maxlag for k in forced_path
        ):
            raise ValueError("forced_path must contain valid P-search lag indices")
        count = len(forced_path)
    else:
        count = int(itmax)
    u0, w0 = np.zeros(nfft), np.zeros(nfft)
    u0[:nt], w0[:nt] = u, w
    du0, dw0 = np.zeros((nfft, npart)), np.zeros((nfft, npart))
    du0[:nt], dw0[:nt] = du, dw
    gaussian = _gauss_filter(dt, nfft, f0)
    uflt, wflt = _gfilter(u0, nfft, gaussian, dt), _gfilter(w0, nfft, gaussian, dt)
    power_u, power_w = np.sum(uflt**2), np.sum(wflt**2)
    if power_u <= 0 or power_w <= 0 or not np.isfinite(power_u + power_w):
        raise ValueError("Filtered waveforms must have nonzero finite energy")
    wf = fft(w0, nfft)
    spikes = np.zeros(nfft)
    dspikes = np.zeros((nfft, npart))

    # For q=dt*spikes, c=a-convolve(q,b), a=C(U,W)/E, b=C(W,W)/E.
    # Precompute da, db once. Each selected c_l updates dc by only two
    # shifted arrays: dc_new=dc-dc_l*roll(b,l)-c_l*roll(db,l).
    if npart:
        duflt = _columns_filter(du0, gaussian, dt)
        dwflt = _columns_filter(dw0, gaussian, dt)
        ufft, wfft = fft(uflt), fft(wflt)
        dufft, dwfft = fft(duflt, axis=0), fft(dwflt, axis=0)
        denergy = 2 * np.sum(wflt[:, None] * dwflt, axis=0)
        a = _correl(uflt, wflt, nfft) / power_w
        b = _correl(wflt, wflt, nfft) / power_w
        da = ifft(
            dufft * wfft.conjugate()[:, None] + ufft[:, None] * dwfft.conjugate(),
            axis=0,
        ).real
        db = ifft(
            dwfft * wfft.conjugate()[:, None] + wfft[:, None] * dwfft.conjugate(),
            axis=0,
        ).real
        dc = (da - a[:, None] * denergy) / power_w
        db = (db - b[:, None] * denergy) / power_w

    residual = uflt
    previous_rms = 1.0
    path, rms, improvements, margins, amplitudes, stop_margins = [], [], [], [], [], []
    degenerate = []
    initial_peak = None
    for iteration in range(count):
        # Primal operation order is deliberately identical to original deconit.
        c = _correl(residual, wflt, nfft) / np.sum(wflt**2)
        candidate = np.abs(c[:maxlag])
        selected = int(np.argmax(candidate)) if forced_path is None else int(forced_path[iteration])
        largest = float(np.max(candidate))
        if initial_peak is None:
            initial_peak = max(largest, np.finfo(float).tiny)
        top_two = np.partition(candidate, -2)[-2:]
        margins.append(float(top_two.max() - top_two.min()))
        amplitude = c[selected]
        path.append(selected)
        amplitudes.append(float(amplitude))
        spikes[selected] = spikes[selected] + amplitude / dt
        if npart:
            damp = dc[selected].copy()  # copy before modifying the selected row
            dspikes[selected] += damp / dt
            # Zero primal residual with a nonzero tangent is a genuine warning;
            # structural Z/Z ties have both primal and tangent near zero.
            if largest < 1e-12 * initial_peak and np.max(np.abs(dc[:maxlag])) > 1e-9:
                degenerate.append(iteration)
            dc -= np.roll(b, selected)[:, None] * damp
            dc -= amplitude * np.roll(db, selected, axis=0)
        predicted = _gfilter(spikes, nfft, gaussian, dt)
        predicted = _gfilter(predicted, nfft, wf, dt)
        residual = uflt - predicted
        current_rms = np.sum(residual**2) / power_u
        improvement = 100 * (previous_rms - current_rms)
        rms.append(float(current_rms))
        improvements.append(float(improvement))
        stop_margins.append(float(abs(improvement) - minderr))
        previous_rms = current_rms
        if forced_path is None and abs(improvement) < minderr:
            break

    rf = _gfilter(spikes, nfft, gaussian, dt)
    shift_samples = int(tshift / dt)  # SeisPy truncates, rather than rounds.
    correction = np.cos(2 * np.pi * shift_samples / nfft)
    if abs(correction) < 1e-8:
        raise ValueError("Singular legacy SeisPy phase-shift normalization")
    rf = _phaseshift(rf, nfft, dt, tshift)[:nt]
    if npart:
        jacobian = _columns_filter(dspikes, gaussian, dt)
        phase = 2 * np.pi * np.arange(1, nfft + 1) * shift_samples / nfft
        phase_factor = np.cos(phase) - 1j * np.sin(phase)
        jacobian = (ifft(fft(jacobian, axis=0) * phase_factor[:, None], axis=0)
                    / correction).real[:nt]
    else:
        jacobian = np.empty((nt, 0))
    if not np.isfinite(rf).all() or not np.isfinite(jacobian).all():
        raise ValueError("Nonfinite iterative result or derivative")
    return IterDeconResult(
        rf, jacobian, path, np.asarray(rms), np.asarray(improvements), np.asarray(margins),
        np.asarray(amplitudes), np.asarray(stop_margins), degenerate, forced_path is not None,
    )


def forward_and_jacobian_iter(
    vp, vs, rho, thickness, *, p=0.06, dt=0.05, npts=1024, shift=10.0, f0=2.0,
    itmax=400, minderr=0.001, periods=None, vp_vs_derivative=None, rho_vs_derivative=None,
    pre_filt=None, zero_halfspace=False,
):
    """Synthetic iterative P-RFs and local fixed-branch Vsapp Jacobian.

    Output columns are independent layer Vs values, including the half-space.
    Optional diagonal Vp/rho slopes give constrained total derivatives.
    Layers/p/filter/iteration tolerances are fixed. Optional pre_filt uses
    the same second-order zero-phase bandpass as SynSeis.filter.
    Diagnostic path stability MUST be checked when interpreting finite changes.

    :param zero_halfspace: Zero the half-space column of all returned Vs
        derivatives without changing the forward response, defaults to False
    :type zero_halfspace: bool, optional
    """
    if not isinstance(zero_halfspace, (bool, np.bool_)):
        raise TypeError('zero_halfspace must be a bool')
    vs = _array(vs, "vs")
    nl = len(vs)
    vp, rho, h = (_array(a, name, nl) for a, name in
                  [(vp, "vp"), (rho, "rho"), (thickness, "thickness")])
    if np.any(vs <= 0) or np.any(vp <= vs) or np.any(rho <= 0):
        raise ValueError("Require positive Vs/rho and Vp > Vs")
    if np.any(h[:-1] <= 0) or h[-1] != 0:
        raise ValueError("Finite layers need positive thickness; half-space thickness must be zero")
    for value, name in [(p, "p"), (dt, "dt"), (f0, "f0"), (shift, "shift")]:
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if np.any(vp * p >= 1):
        raise ValueError("Only subcritical incidence p*Vp < 1 is supported")
    if isinstance(npts, (bool, np.bool_)) or not isinstance(npts, (int, np.integer)) or npts < 8:
        raise ValueError("npts must be an integer >= 8")
    if f0 * dt > 0.5:
        raise ValueError("Gaussian pulse is undersampled")
    periods = _array(np.linspace(0.1, 3.5, 35) if periods is None else periods, "periods")
    if np.any(periods <= 0) or np.any(np.diff(periods) <= 0) or periods[0] < 2 * dt - 1e-9:
        raise ValueError("periods must increase strictly and be >= 2*dt")
    nfft = 1 << (int(npts) - 1).bit_length()
    times = np.arange(nfft) * dt - shift
    if times[0] > -periods[-1] + 1e-8 or times[-1] < periods[-1] - 1e-8:
        raise ValueError("Waveform span must cover every [-T,T] window")
    da = np.zeros(nl) if vp_vs_derivative is None else _array(vp_vs_derivative, "dVp/dVs", nl)
    dr = np.zeros(nl) if rho_vs_derivative is None else _array(rho_vs_derivative, "dRho/dVs", nl)
    ur, uz, dur, duz = _spectra_and_tangents(vp, vs, rho, h, da, dr, p, dt, nfft)
    r, z = ifft(ur).real[::-1] / nfft, -ifft(uz).real[::-1] / nfft
    d_r, d_z = ifft(dur, axis=0).real[::-1] / nfft, -ifft(duz, axis=0).real[::-1] / nfft
    r, z, d_r, d_z = _prefilter_waveforms(r, z, d_r, d_z, dt, pre_filt)
    options = {"tshift": shift, "f0": f0, "itmax": itmax, "minderr": minderr}
    radial = deconit_tangent(r, z, d_r, d_z, dt, **options)
    vertical = deconit_tangent(z, z, d_z, d_z, dt, **options)
    curve = _compute_vsapp_from_components(times, radial.rf, vertical.rf, rayp=p, periods_s=periods)
    if any(status != "ok" for status in curve.status):
        raise ValueError(f"Vsapp is outside the valid measurement branch: {curve.status}")
    operator = _integration_operator(times, periods)
    a, b = curve.radial_integral, curve.vertical_integral
    da_int = np.einsum("it,tj->ij", operator, radial.jacobian, optimize=False)
    db_int = np.einsum("it,tj->ij", operator, vertical.jacobian, optimize=False)
    factor = np.cos(0.5 * np.arctan(a / b)) / (2 * p * (a * a + b * b))
    jacobian = factor[:, None] * (b[:, None] * da_int - a[:, None] * db_int)
    diagnostics = {
        "deconvolution": "iter", "derivative_method": "analytic, fixed execution branch",
        "nfft": nfft, "n_layers_including_halfspace": nl,
        "prefilter": None if pre_filt is None else list(pre_filt),
        "itmax": int(itmax), "minderr": float(minderr),
        "actual_phase_shift_s": int(shift / dt) * dt,
        "branch_warning": bool(radial.degenerate_tangent_iterations
                               or vertical.degenerate_tangent_iterations),
    }
    for name, result in [("radial", radial), ("vertical", vertical)]:
        diagnostics.update({
            f"{name}_path": result.path,
            f"{name}_iterations": len(result.path),
            f"{name}_rms": result.rms.tolist(),
            f"{name}_improvements": result.improvements.tolist(),
            f"{name}_peak_margins": result.peak_margins.tolist(),
            f"{name}_amplitudes": result.amplitudes.tolist(),
            f"{name}_stop_margins": result.stop_margins.tolist(),
            f"{name}_degenerate_tangent_iterations": result.degenerate_tangent_iterations,
        })
    diagnostics['zero_halfspace'] = bool(zero_halfspace)
    if zero_halfspace:
        jacobian[:, -1] = 0.
        radial.jacobian[:, -1] = 0.
        vertical.jacobian[:, -1] = 0.
    return VsappKernelResult(times, radial.rf, vertical.rf, curve.vs_km_s,
                             radial.jacobian, vertical.jacobian, jacobian, diagnostics,
                             periods.copy(), h.copy())

def _prefilter_waveforms(r, z, dr, dz, dt, pre_filt):
    """Apply the fixed SynSeis.filter bandpass to primal and tangent arrays."""
    if pre_filt is None:
        return r, z, dr, dz
    from obspy.signal.filter import bandpass

    limits = _array(pre_filt, "pre_filt", 2)
    if not 0 < limits[0] < limits[1] < 0.5 / dt:
        raise ValueError("pre_filt must satisfy 0 < freqmin < freqmax < Nyquist")

    def apply(trace):
        return bandpass(trace, limits[0], limits[1], df=1.0 / dt,
                        corners=2, zerophase=True)

    # Filter 1-D columns to preserve older ObsPy's zero-phase axis semantics.
    return (apply(r), apply(z),
            np.column_stack([apply(column) for column in dr.T]),
            np.column_stack([apply(column) for column in dz.T]))


def vsapp_kernel(depmod, rayp, periods_s, *, method="iter", zero_halfspace=False, **kwargs):
    """Return RFs, Vsapp and all Vs derivatives for an existing DepModel.

    rayp is in s/km; periods_s are cosine-squared window half-widths in s.
    method is 'iter' (default) or 'water'. Other keywords are forwarded to
    the matching array-level function (dt, npts, shift, f0, pre_filt and
    method-specific settings). pre_filt=None means no front-end bandpass.
    Kernel columns correspond to depmod.vs, including any grid resampling
    performed by DepModel.read_layer_model. The model is never modified.

    :param zero_halfspace: Set the last columns of jacobian, dradial_dvs and
        dvertical_dvs to zero. Forward responses are unchanged, defaults to False
    :type zero_halfspace: bool, optional
    :return: Synthetic RFs, apparent Vs and a plottable sensitivity matrix
    :rtype: seispy.vsapp_kernel.VsappKernelResult
    """
    from seispy.core.depmodel import DepModel

    if not isinstance(depmod, DepModel):
        raise TypeError("depmod must be a seispy.DepModel")
    if method not in ("iter", "water"):
        raise ValueError("method must be 'iter' or 'water'")
    function = forward_and_jacobian_iter if method == "iter" else forward_and_jacobian
    result = function(depmod.vp, depmod.vs, depmod.rho, depmod.thickness,
                      p=rayp, periods=periods_s, zero_halfspace=zero_halfspace, **kwargs)
    result.diagnostics.update({
        "rayp_s_km": float(rayp),
        "periods_s": np.asarray(periods_s, dtype=float).tolist(),
        "model_parameters": "current DepModel.vs entries, including half-space",
    })
    return result
