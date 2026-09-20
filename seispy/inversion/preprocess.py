# SPDX-License-Identifier: GPL-3.0-only
"""Validate inversion inputs, assemble constraints and select station events."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass

import numpy as np

from seispy.core.depmodel import DepModel
from seispy.vsapp import StationVsAppResult, period_grid, positive, vector

from .gradient import _uncertainty


@dataclass(frozen=True)
class InversionProblem:
    """Owned data, model template and fixed settings for one inversion."""

    template: DepModel
    periods: np.ndarray
    observed: np.ndarray
    uncertainty: np.ndarray
    rayp: float
    f0: float
    method: str
    options: dict
    initial_vs: np.ndarray
    depths_km: np.ndarray
    fixed: np.ndarray
    zero_halfspace: bool


def prepare_inversion(initial_model, periods_s, observed_vsapp, *, rayp, f0, sigma,
                      method, zero_halfspace, kernel_kwargs):
    """Create an independent problem without executing a forward calculation."""
    if not isinstance(initial_model, DepModel):
        raise TypeError('initial_model must be a DepModel')
    # build periods and observed data
    periods = period_grid(periods_s)
    observed = vector(observed_vsapp, 'observed_vsapp')
    if observed.shape != periods.shape or np.any(observed <= 0):
        raise ValueError('observed_vsapp must be positive and match periods_s')
    uncertainty = _uncertainty(sigma, observed.size)
    rayp, f0 = positive(rayp, 'rayp'), positive(f0, 'f0')
    if method not in ('iter', 'water'):
        raise ValueError("method must be 'iter' or 'water'")
    if not isinstance(zero_halfspace, (bool, np.bool_)):
        raise TypeError('zero_halfspace must be a bool')
    options = _kernel_options(kernel_kwargs, method)
    template = deepcopy(initial_model)
    initial_vs = vector(template.vs, 'initial_model.vs')
    thickness = vector(template.thickness, 'initial_model.thickness')
    if (thickness.shape != initial_vs.shape or np.any(thickness[:-1] <= 0)
            or thickness[-1] != 0):
        raise ValueError('Model needs positive finite-layer thicknesses and a final half-space')
    depths = np.r_[0., np.cumsum(thickness[:-1])] + thickness / 2
    if np.any(initial_vs <= 0):
        raise ValueError('Initial Vs must be positive for ln(Vs) optimization')
    # zero_halfspace masks only the kernel column. Keep the half-space free
    # so depth smoothing can propagate neighbouring updates into it.
    fixed = np.zeros(initial_vs.size, dtype=bool)
    return InversionProblem(template, periods, observed, uncertainty, rayp, f0, method,
                             options, initial_vs, depths, fixed,
                             bool(zero_halfspace))


def prepare_station_data(measurements):
    """Return periods, mean curve, mean slowness and a whole-event QC mask."""
    if not isinstance(measurements, StationVsAppResult):
        raise TypeError('measurements must be a StationVsAppResult')
    values = np.asarray(measurements.vs_km_s, dtype=float)
    status = np.asarray(measurements.status)
    rayps = vector(measurements.rayp, 'measurements.rayp')
    periods = period_grid(measurements.periods_s)
    if (values.shape != (rayps.size, periods.size) or status.shape != values.shape
            or np.any(rayps <= 0)):
        raise ValueError('Station curves, statuses, periods and ray parameters must match')
    valid = np.all((status == 'ok') & np.isfinite(values) & (values > 0), axis=1)
    if not np.any(valid):
        raise ValueError('No event is valid at every requested period')
    return periods, values[valid].mean(axis=0), float(rayps[valid].mean()), valid


def _kernel_options(kernel_kwargs, method):
    """Own fixed processing options and disallow overriding material coupling."""
    if kernel_kwargs is not None and not isinstance(kernel_kwargs, Mapping):
        raise TypeError('kernel_kwargs must be a mapping')
    options = deepcopy(dict(kernel_kwargs or {}))
    allowed = {'dt', 'npts', 'shift', 'pre_filt'}
    allowed |= {'itmax', 'minderr'} if method == 'iter' else {'wlevel'}
    if options.keys() - allowed:
        raise TypeError('Unsupported kernel_kwargs: '
                        + ', '.join(sorted(map(str, options.keys() - allowed))))
    return options
