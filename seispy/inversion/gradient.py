# SPDX-License-Identifier: GPL-3.0-only
"""Weighted Vsapp data gradients and Gaussian smoothing in physical depth."""

import numpy as np

from seispy.signal import smooth
from seispy.vsapp import vector
from seispy.vsapp_kernel import VsappKernelResult


def vsapp_gradient(kernel, observed_vsapp, *, vs_km_s, sigma=None):
    """Compute a weighted least-squares misfit and its unsmoothed ln(Vs) gradient.

    For N periods, misfit = sum(((predicted - observed) / sigma)**2) / (2*N)
    and gradient = Vs * (J.T @ ((predicted - observed) / sigma**2)) / N,
    where J is the kernel derivative with respect to Vs. Thus each kernel
    column is multiplied by the current layer Vs before depth smoothing. No
    normalization of kernel columns is applied. Zeroed half-space columns
    in the input kernel remain zero in the resulting gradient.

    :param kernel: Forward response and layer Jacobian
    :type kernel: seispy.vsapp_kernel.VsappKernelResult
    :param vs_km_s: Current positive layer Vs in km/s, matching the kernel model
    :type vs_km_s: array_like
    :param observed_vsapp: Positive measured apparent velocities in km/s
    :type observed_vsapp: array_like
    :param sigma: Positive data uncertainties in km/s, scalar or one per period.
        None gives equal unit weights, defaults to None
    :type sigma: float or array_like, optional
    :return: Misfit and one derivative with respect to ln(Vs) per model layer
    :rtype: (float, numpy.ndarray)
    :raises ValueError: If data, uncertainties or kernel dimensions are invalid
    :raises TypeError: If kernel is not a VsappKernelResult
    """
    if not isinstance(kernel, VsappKernelResult):
        raise TypeError('kernel must be a VsappKernelResult')
    observed = vector(observed_vsapp, 'observed_vsapp')
    predicted = vector(kernel.vsapp, 'kernel.vsapp')
    if predicted.shape != observed.shape or np.any(observed <= 0):
        raise ValueError('observed_vsapp must be positive and match kernel.vsapp')
    uncertainty = _uncertainty(sigma, observed.size)
    jacobian = np.asarray(kernel.jacobian, dtype=float)
    if (jacobian.ndim != 2 or jacobian.shape[0] != observed.size
            or jacobian.shape[1] == 0 or not np.isfinite(jacobian).all()):
        raise ValueError('kernel.jacobian must be finite with shape (periods, layers)')
    speed = vector(vs_km_s, 'vs_km_s')
    if speed.shape != (jacobian.shape[1],) or np.any(speed <= 0):
        raise ValueError('vs_km_s must be positive and match kernel layers')
    residual = (predicted - observed) / uncertainty
    misfit = float(np.sum(residual**2) / (2 * observed.size))
    gradient = speed * np.einsum('ij,i->j', jacobian, residual / uncertainty) / observed.size
    if not np.isfinite(misfit) or not np.isfinite(gradient).all():
        raise ValueError('Nonfinite weighted misfit or gradient; check data uncertainties')
    return misfit, gradient


def smooth_gradient(gradient, depths_km, sigma_km, *, fixed_mask=None):
    """Gaussian-smooth a model gradient with :func:`seispy.signal.smooth`.

    Convert sigma_km to samples and truncate the Gaussian at four standard
    deviations. Uneven grids are linearly resampled at a spacing no larger
    than their smallest depth interval, then interpolated back. Boundaries
    are reflected. Smooth the free-parameter mask as well to exclude fixed
    entries without diluting nearby gradients. Fixed entries remain zero.
    A zero width disables smoothing. Inputs are not modified.

    This normalized smoothing is a search-direction preconditioner,
    not a model smoothness penalty or the exact objective gradient. The
    L-BFGS checks descent and may retry without smoothing. SD applies its
    gradient-angle-controlled step directly, with no line search.

    :param gradient: One finite gradient value per model layer
    :type gradient: array_like
    :param depths_km: Strictly increasing layer representative depths in km
    :type depths_km: array_like
    :param sigma_km: Nonnegative Gaussian standard deviation in km
    :type sigma_km: float
    :param fixed_mask: Boolean mask of excluded model parameters, defaults to None
    :type fixed_mask: array_like, optional
    :return: Smoothed gradient with fixed entries set to zero
    :rtype: numpy.ndarray
    :raises ValueError: If coordinates, width or fixed mask are invalid
    """
    values = vector(gradient, 'gradient')
    depths = vector(depths_km, 'depths_km')
    if depths.shape != values.shape or np.any(np.diff(depths) <= 0):
        raise ValueError('depths_km must increase strictly and match gradient')
    width = float(sigma_km)
    if not np.isfinite(width) or width < 0:
        raise ValueError('sigma_km must be finite and nonnegative')
    fixed = np.zeros(values.size, dtype=bool)
    if fixed_mask is not None:
        fixed = np.asarray(fixed_mask)
        if fixed.dtype.kind != 'b' or fixed.shape != values.shape:
            raise ValueError('fixed_mask must be boolean and match gradient')
    values[fixed] = 0.
    indices = np.flatnonzero(~fixed)
    if width == 0 or indices.size < 2:
        return values
    # Convolution needs a uniform depth grid; the public width stays in km.
    spacing = np.diff(depths)
    grid = depths
    if not np.allclose(spacing, spacing[0], rtol=1e-6, atol=0.):
        count = int(np.ceil((depths[-1] - depths[0]) / spacing.min())) + 1
        grid = np.linspace(depths[0], depths[-1], count)
    sigma_samples = width / (grid[1] - grid[0])
    options = dict(half_len=int(np.ceil(4 * sigma_samples)),
                   window='gaussian', sigma=sigma_samples)
    numerator = smooth(np.interp(grid, depths, values), **options)
    denominator = smooth(np.interp(grid, depths, (~fixed).astype(float)), **options)
    filtered = np.divide(numerator, denominator, out=np.zeros_like(numerator),
                         where=denominator > 0)
    result = np.interp(depths, grid, filtered)
    result[fixed] = 0.
    return result


def _uncertainty(sigma, size):
    """Validate scalar or per-period standard deviations."""
    values = np.ones(size) if sigma is None else np.asarray(sigma, dtype=float)
    if values.ndim == 0:
        values = np.full(size, float(values))
    if values.shape != (size,) or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('sigma must be finite, positive and scalar or one value per period')
    return values.copy()
