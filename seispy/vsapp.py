# SPDX-License-Identifier: GPL-3.0-only
"""Array-only implementation of Yao et al. (2022), equations (2)-(3).

Time is in seconds relative to the direct P arrival, p in s/km and Vs in
km/s. Public measurements take radial RFs and construct the reference internally.
"""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray


class VsAppError(ValueError):
    """Invalid waveform, measurement window, or Vsapp parameter."""

FloatArray = NDArray[np.float64]


def vector(value: ArrayLike, name: str, *, minimum_size: int = 1) -> FloatArray:
    """Make an owned, finite, one-dimensional float array."""
    try:
        array = np.array(value, dtype=np.float64, copy=True)
    except (ValueError, TypeError) as exc:
        raise VsAppError(f"{name} must contain real numbers") from exc
    if array.ndim != 1 or array.size < minimum_size or not np.isfinite(array).all():
        raise VsAppError(f"{name} must be finite, 1-D, and have >= {minimum_size} samples")
    return array


def positive(value: float, name: str) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise VsAppError(f"{name} must be a finite positive number") from exc
    if not np.isfinite(value) or value <= 0:
        raise VsAppError(f"{name} must be finite and > 0")
    return value


def time_axis(time_s: ArrayLike) -> FloatArray:
    times = vector(time_s, "time_s", minimum_size=3)
    steps = np.diff(times)
    if np.any(steps <= 0) or not np.allclose(steps, steps[0], rtol=1e-5, atol=1e-9):
        raise VsAppError("time_s must be strictly increasing and uniformly sampled")
    if times[0] > 0 or times[-1] < 0:
        raise VsAppError("time_s must include the direct P arrival (t=0)")
    return times


def period_grid(periods_s: ArrayLike) -> FloatArray:
    periods = vector(periods_s, "periods_s")
    if np.any(periods <= 0) or np.any(np.diff(periods) <= 0):
        raise VsAppError("periods_s must be positive and strictly increasing; Vs0 is separate")
    return periods


@dataclass(frozen=True)
class Curve:
    """Measurements and QC; rayp is in s/km. Invalid velocities remain NaN."""

    periods_s: FloatArray
    radial_integral: FloatArray
    vertical_integral: FloatArray
    ratio: FloatArray
    angle_deg: FloatArray
    vs_km_s: FloatArray
    status: tuple[str, ...]
    radial_zero: float
    vertical_zero: float
    vs0_km_s: float
    vs0_status: str
    rayp: float


def gaussian_vertical(
    time_s: ArrayLike,
    gaussian_factor: float,
    *,
    normalization: Literal["area", "peak"],
) -> FloatArray:
    """Ideal Z/Z reference for G(omega)=exp[-omega**2/(4*a**2)].

    area: a/sqrt(pi)*exp[-(a*t)**2], integral 1 (unnormalized deconvolution).
    peak: exp[-(a*t)**2], peak 1 (R must have been divided by the same peak).
    This is an idealization, not an observed vertical RF or water-level correction.
    """
    times = time_axis(time_s)
    factor = positive(gaussian_factor, "gaussian_factor")
    if normalization not in ("area", "peak"):
        raise VsAppError("normalization must be 'area' or 'peak'")
    if factor * np.median(np.diff(times)) > 0.5:
        raise VsAppError("Gaussian pulse is undersampled (a*dt > 0.5)")
    pulse = np.exp(-np.square(factor * times))
    return pulse * (factor / np.sqrt(np.pi) if normalization == "area" else 1.0)


def seispy_iter_vertical(
    time_s: ArrayLike, gaussian_factor: float, *, deconvolution_npts: int | None = None,
    deconvolution_shift_samples: int | None = None,
) -> FloatArray:
    """Reproduce the discrete SeisPy iterative Z/Z impulse response.

    Matches decon.gaussFilter + phaseshift, including its legacy phase-shift
    amplitude factor. Valid only for unresampled, unscaled iterative RFs from
    that implementation. npts defaults to the saved trace length. If the
    deconvolution input had another length, supply it explicitly.

    Without ``deconvolution_shift_samples``, reconstruct the nominal P-aligned
    sample. To reproduce native ``phaseshift`` at a floating-point boundary,
    supply ``int(original_tshift / original_dt)`` instead. For example,
    ``int(0.3 / 0.1)`` is 2, although the nominal aligned sample is 3.
    SAC headers round their timing values; do not infer the original index
    with ``int`` from SAC's approximate delta. An explicit index also allows
    a non-grid-aligned nominal P arrival. The RF and reference retain the
    original discrete shift; no additional alignment is applied.
    """
    times = time_axis(time_s)
    factor = positive(gaussian_factor, "gaussian_factor")
    dt = float(np.median(np.diff(times)))
    nt = len(times) if deconvolution_npts is None else deconvolution_npts
    if isinstance(nt, bool) or not isinstance(nt, (int, np.integer)) or nt < len(times):
        raise VsAppError("deconvolution_npts must be an integer >= saved trace length")
    nfft = 1 << (int(nt) - 1).bit_length()
    if deconvolution_shift_samples is None:
        shift_samples = -times[0] / dt
        shift = int(round(shift_samples))
        if not np.isclose(shift_samples, shift, atol=1e-5, rtol=0):
            raise VsAppError('SeisPy iterative reference requires P aligned to a sample '
                             'or an explicit deconvolution_shift_samples')
    else:
        if (isinstance(deconvolution_shift_samples, (bool, np.bool_))
                or not isinstance(deconvolution_shift_samples, (int, np.integer))
                or not 0 <= deconvolution_shift_samples < nfft):
            raise VsAppError('deconvolution_shift_samples must be an integer in [0, nfft)')
        shift = int(deconvolution_shift_samples)
    gaussian_vertical(times, factor, normalization="area")  # resolution validation
    omega = 2 * np.pi * np.fft.fftfreq(nfft, dt)
    spectrum = np.exp(-0.25 * (omega / factor) ** 2) / dt
    phase = 2 * np.pi * np.arange(1, nfft + 1) * shift / nfft
    correction = np.cos(2 * np.pi * shift / nfft)
    if abs(correction) < 1e-8:
        raise VsAppError("Singular legacy SeisPy phase-shift correction; "
                         "recompute the RF with a different deconvolution shift")
    return (np.fft.ifft(spectrum * np.exp(-1j * phase)).real / correction)[: len(times)]


def _velocity(
    radial: float, vertical: float, floor: float, p: float
) -> tuple[float, float, float, str]:
    if not np.isfinite(radial) or not np.isfinite(vertical):
        return np.nan, np.nan, np.nan, "nonfinite_integral"
    if vertical <= floor:
        return np.nan, np.nan, np.nan, "invalid_vertical"
    ratio = radial / vertical
    if not np.isfinite(ratio):
        return np.nan, np.nan, np.nan, "nonfinite_ratio"
    angle = np.arctan(ratio)
    if radial <= 0:
        return ratio, np.degrees(angle), np.nan, "nonpositive_radial"
    with np.errstate(over="ignore"):
        speed = np.sin(angle / 2) / p
    if not np.isfinite(speed):
        return ratio, np.degrees(angle), np.nan, "nonfinite_velocity"
    return ratio, np.degrees(angle), speed, "ok"


def compute_vsapp(
    time_s: ArrayLike,
    radial: ArrayLike,
    *,
    rayp: float,
    periods_s: ArrayLike,
    f0: float,
    reference: str = "seispy-iter",
    deconvolution_npts: int | None = None,
    deconvolution_shift_samples: int | None = None,
    denominator_rtol: float = 1e-10,
) -> Curve:
    """Measure Vs0 and Vsapp(T) from a radial P receiver function.

    ``rayp`` is in s/km; ``time_s`` is relative to direct P and ``periods_s``
    contains cosine-squared half-window widths in seconds. ``f0`` is the
    Gaussian factor used during deconvolution. No Z/Z waveform is required.

    The default ``reference='seispy-iter'`` constructs the discrete reference
    for original SeisPy iterative RFs with no subsequent resampling, filtering
    or independent amplitude normalization. For tail-truncated RFs, supply
    the original input length as ``deconvolution_npts``. An explicit
    ``deconvolution_shift_samples=int(original_tshift / original_dt)`` retains
    SeisPy's integer truncation at floating-point sample boundaries; otherwise
    the nominal P-aligned sample is used. Do not infer that index by truncating
    rounded SAC timing headers.

    ``'gaussian-area'`` and ``'gaussian-peak'`` select ideal unit-area and
    unit-peak references, respectively. For the latter, R must have been
    divided by the same vertical peak. These are declared amplitude
    conventions, not corrections for water-level deconvolution. Input arrays
    are not modified. See ``Curve`` for integrals and QC alongside velocities.
    """
    if reference not in ('seispy-iter', 'gaussian-area', 'gaussian-peak'):
        raise VsAppError("reference must be 'seispy-iter', 'gaussian-area', or 'gaussian-peak'")
    if (reference != 'seispy-iter'
            and (deconvolution_npts is not None or deconvolution_shift_samples is not None)):
        raise VsAppError('deconvolution_npts and deconvolution_shift_samples '
                         'are only valid with reference=seispy-iter')
    factor = positive(f0, 'f0')
    if reference == 'seispy-iter':
        vertical = seispy_iter_vertical(
            time_s, factor, deconvolution_npts=deconvolution_npts,
            deconvolution_shift_samples=deconvolution_shift_samples,
        )
    else:
        vertical = gaussian_vertical(
            time_s, factor, normalization=reference.removeprefix('gaussian-'),
        )
    return _compute_vsapp_from_components(
        time_s, radial, vertical, rayp=rayp, periods_s=periods_s,
        denominator_rtol=denominator_rtol,
    )


def _compute_vsapp_from_components(
    time_s: ArrayLike,
    radial: ArrayLike,
    vertical: ArrayLike,
    *,
    rayp: float,
    periods_s: ArrayLike,
    denominator_rtol: float = 1e-10,
) -> Curve:
    """Internal paired-component integration for synthetic RFs and kernels.

    R and Z must share a deconvolution, filter and amplitude scale.
    T is the half-width of a cosine-squared kernel.
    rayp is the event's horizontal slowness in s/km. time_s is relative to its
    direct P arrival; this function does not align peaks or apply moveout.
    Both components are smoothed on [-T,T], then their zero-lag amplitude
    ratio is converted with sin(arctan(R/Z)/2)/p. Numerical integration uses
    trapezoids at original samples plus interpolated endpoints and t=0.
    Requests outside the data window or T < 2*dt are rejected, not padded.
    Negative R integrals are retained in diagnostics but produce NaN Vs.
    """
    times = time_axis(time_s)
    r = vector(radial, "radial", minimum_size=3)
    z = vector(vertical, "vertical", minimum_size=3)
    if r.shape != times.shape or z.shape != times.shape:
        raise VsAppError("time_s, radial and vertical must have identical shapes")
    p = positive(rayp, "rayp")
    tolerance = positive(denominator_rtol, "denominator_rtol")
    if tolerance >= 1:
        raise VsAppError("denominator_rtol must be < 1")
    periods = period_grid(periods_s)
    dt = float(np.median(np.diff(times)))
    if periods[0] < 2 * dt - 1e-9:
        raise VsAppError(f"Minimum T must be >= 2*dt ({2 * dt:g} s)")
    if times[0] > -periods[-1] + 1e-8 or times[-1] < periods[-1] - 1e-8:
        raise VsAppError(f"RF must cover [-{periods[-1]:g}, {periods[-1]:g}] s")
    z_scale = float(np.max(np.abs(z)))
    if z_scale == 0:
        raise VsAppError("Vertical RF is identically zero")
    r0, z0 = float(np.interp(0, times, r)), float(np.interp(0, times, z))
    _, _, vs0, status0 = _velocity(r0, z0, z_scale * tolerance, p)
    r_int, z_int, ratios, angles, speeds, statuses = [], [], [], [], [], []
    for period in periods:
        inside = times[(times > -period) & (times < period)]
        grid = np.unique(np.concatenate(([-period, 0.0, period], inside)))
        weight = np.cos(np.pi * grid / (2 * period)) ** 2
        weight[[0, -1]] = 0.0
        # Explicit trapezoid sum supports NumPy 1.24 and 2.x without deprecated APIs.
        dr = np.interp(grid, times, r) * weight
        dz = np.interp(grid, times, z) * weight
        widths = np.diff(grid)
        ri = float(np.sum((dr[1:] + dr[:-1]) * widths / 2))
        zi = float(np.sum((dz[1:] + dz[:-1]) * widths / 2))
        ratio, angle, speed, status = _velocity(ri, zi, z_scale * period * tolerance, p)
        r_int.append(ri)
        z_int.append(zi)
        ratios.append(ratio)
        angles.append(angle)
        speeds.append(speed)
        statuses.append(status)
    arrays = [periods, *map(np.asarray, (r_int, z_int, ratios, angles, speeds))]
    for array in arrays:
        array.setflags(write=False)
    return Curve(*arrays, tuple(statuses), r0, z0, vs0, status0, p)


@dataclass(frozen=True)
class StationVsAppResult:
    """Event-by-event Vsapp measurements in the station's current order.

    ``rayp`` is in s/km and velocities are in km/s. ``vs_km_s`` and
    ``status`` have shape (events, periods). Array fields are independent,
    read-only snapshots: later station sorting or processing cannot change
    this result. ``curves`` retains each event's integrals and full QC.
    """

    curves: tuple[Curve, ...]
    event: NDArray[np.str_]
    rayp: FloatArray
    periods_s: FloatArray
    vs_km_s: FloatArray
    vs0_km_s: FloatArray
    status: NDArray[np.str_]
    vs0_status: tuple[str, ...]
    reference: str


def compute_station_vsapp(
    station,
    periods_s: ArrayLike,
    *,
    reference: str = "seispy-iter",
    deconvolution_npts: int | None = None,
    deconvolution_shift_samples: int | None = None,
    denominator_rtol: float = 1e-10,
) -> StationVsAppResult:
    """Measure radial P receiver functions held by a SeisPy RFStation.

    Uses only ``station.data_prime`` and existing metadata. Each event's
    ``f0`` sets the internally constructed reference:

    * ``'seispy-iter'`` (default) reconstructs SeisPy's discrete iterative response.
      It requires native, unresampled RFs whose radial amplitudes have not
      been independently normalized or filtered. If the saved RF was tail-truncated, supply
      the original deconvolution input length as ``deconvolution_npts``.
      That length and prior resampling/scaling cannot be inferred reliably
      from an RFStation or SAC headers.
    * ``'gaussian-area'`` assumes an ideal unit-area Gaussian Z/Z pulse.
    * ``'gaussian-peak'`` assumes an ideal unit-peak Gaussian and requires R
      to have been divided by the same vertical peak, not its own peak.

    The Gaussian choices are idealizations, not a water-level correction.
    No vertical waveform input is needed. No normalization,
    moveout correction, stacking or station mutation is performed here.
    Each event uses its own slowness (converted from RFStation's s/rad to
    s/km). Averaging slownesses or stacking waveforms before conversion
    generally produces a different measurement.

    Native iterative RFTrace streams preserve their discrete shift while
    station sampling and shift are unchanged. SAC input has no such timing
    provenance and defaults to nominal P alignment. If the original
    ``int(tshift / dt)`` differs from the nominal aligned sample, supply it
    as ``deconvolution_shift_samples``. Use the
    original deconvolution timing, not rounded SAC headers, for that index.
    """
    from seispy.geo import srad2skm

    if station.comp != 'R':
        raise VsAppError("Vsapp requires radial P receiver functions (comp='R')")
    if getattr(station, 'prime_phase', '').upper() == 'S':
        raise VsAppError('Vsapp requires P, not S, receiver functions')
    times = time_axis(station.time_axis)
    try:
        radial = np.array(station.data_prime, dtype=np.float64, copy=True)
    except (TypeError, ValueError) as exc:
        raise VsAppError('station.data_prime must contain real numbers') from exc
    if (radial.ndim != 2 or radial.shape[0] == 0
            or radial.shape[1] != times.size or not np.isfinite(radial).all()):
        raise VsAppError('station.data_prime must be finite with shape (events, samples)')
    n_events = radial.shape[0]
    if station.ev_num != n_events:
        raise VsAppError('station.ev_num does not match station.data_prime')
    events = np.array(station.event, dtype=str, copy=True)
    phases = np.asarray(station.phase, dtype=str)
    if events.shape != (n_events,) or phases.shape != (n_events,):
        raise VsAppError('station event and phase metadata must match its waveforms')
    if any(phase.upper().endswith('S') for phase in phases):
        raise VsAppError('Vsapp requires P, not S, receiver functions')
    rayp = srad2skm(vector(station.rayp, 'station.rayp'))
    if rayp.shape != (n_events,) or np.any(rayp <= 0):
        raise VsAppError('station.rayp must contain one positive slowness per event')
    periods = period_grid(periods_s)

    factors = vector(station.f0, 'station.f0')
    if factors.shape != (n_events,) or np.any(factors <= 0):
        raise VsAppError('station.f0 must contain one positive Gaussian factor per event')
    if reference == 'seispy-iter':
        timing = getattr(station, '_vsapp_iter_timing', None)
        if (deconvolution_shift_samples is None and timing is not None
                and station.sampling == timing[0] and station.shift == timing[1]
                and np.array_equal(times, np.arange(times.size) * timing[0] - timing[1])):
            deconvolution_shift_samples = timing[2]
    curves = tuple(
        compute_vsapp(times, r, rayp=p, periods_s=periods, f0=factor,
                      reference=reference, deconvolution_npts=deconvolution_npts,
                      deconvolution_shift_samples=deconvolution_shift_samples,
                      denominator_rtol=denominator_rtol)
        for r, p, factor in zip(radial, rayp, factors, strict=True)
    )
    velocities = np.stack([curve.vs_km_s for curve in curves])
    vs0 = np.asarray([curve.vs0_km_s for curve in curves])
    statuses = np.asarray([curve.status for curve in curves], dtype=str)
    for array in (events, rayp, periods, velocities, vs0, statuses):
        array.setflags(write=False)
    return StationVsAppResult(
        curves, events, rayp, periods, velocities, vs0, statuses,
        tuple(curve.vs0_status for curve in curves), reference,
    )


__all__ = [
    'VsAppError', 'Curve', 'StationVsAppResult', 'compute_vsapp',
    'compute_station_vsapp', 'gaussian_vertical', 'seispy_iter_vertical',
]
