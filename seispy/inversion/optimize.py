# SPDX-License-Identifier: GPL-3.0-only
"""Stopping criteria and accepted-state summary; settings live in InvPara."""

from dataclasses import dataclass

import numpy as np

from seispy.core.depmodel import DepModel

from .forward import ForwardState
from .postprocess import VsappIteration


@dataclass(frozen=True)
class OptimizationOutcome:
    """Final accepted state and trace returned to the postprocessor."""

    state: ForwardState
    initial_model: DepModel
    history: tuple[VsappIteration, ...]
    converged: bool
    message: str
    n_evaluations: int


def gradient_norm(state):
    """Infinity norm of the raw ln(Vs) gradient, with fixed entries zeroed."""
    return float(np.max(np.abs(state.gradient)))


def gradient_angle(previous, current):
    """Return the angle in degrees between two raw ln(Vs) gradients.

    Fixed coordinates have already been zeroed by the forward evaluator.
    Return None for the first iteration or if either gradient is zero.
    Scale before normalization to avoid overflow or underflow in the norm.

    :param previous: Previous accepted model's gradient, or None initially
    :type previous: numpy.ndarray or None
    :param current: Current accepted model's gradient
    :type current: numpy.ndarray
    :return: Angle in [0, 180] degrees, or None when undefined
    :rtype: float or None
    """
    if previous is None:
        return None
    previous_scale = float(np.max(np.abs(previous)))
    current_scale = float(np.max(np.abs(current)))
    if previous_scale == 0 or current_scale == 0:
        return None
    previous_unit = previous / previous_scale
    previous_unit /= np.linalg.norm(previous_unit)
    current_unit = current / current_scale
    current_unit /= np.linalg.norm(current_unit)
    cosine = float(np.sum(previous_unit * current_unit))
    return float(np.degrees(np.arccos(np.clip(cosine, -1., 1.))))


def windowed_misfit_change(history, window):
    """Compare two adjacent means of accepted-iteration misfits.

    Exclude iteration zero. For the last 2*window accepted updates, let B
    be the earlier mean and A the later mean. Return abs(A-B)/B, or None
    until both windows are full. If B=0, return zero when A=0 and infinity
    otherwise. Rejected line-search trials are absent from history.

    :param history: Accepted iteration records, including iteration zero
    :type history: sequence of VsappIteration
    :param window: Number of accepted updates per window; positive integer
    :type window: int
    :return: Relative change of window means, or None if there is insufficient history
    :rtype: float or None
    """
    if len(history) - 1 < 2 * window:
        return None
    misfits = np.array([item.misfit for item in history[-2 * window:]])
    scale = float(np.max(misfits))
    if scale == 0:
        return 0.
    # A common scale cancels in the ratio and keeps the means numerically safe.
    misfits /= scale
    previous_mean = float(np.mean(misfits[:window]))
    recent_mean = float(np.mean(misfits[window:]))
    if previous_mean == 0:
        return float('inf')
    return abs(recent_mean - previous_mean) / previous_mean
