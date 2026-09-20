# SPDX-License-Identifier: GPL-3.0-only
"""Wolfe search with a decreasing-step fallback at the L-BFGS update cap."""

import numpy as np


def strong_wolfe(problem, controls, forward, state, direction):
    """Find a strong Wolfe step, or a decreasing step at the update cap.

    ``direction`` is the proposed descent direction in ln(Vs), including
    any smoothing. All directional derivatives use the raw gradient from
    the kernel, including its optional half-space column masking.
    Start at alpha=min(1, amax), expand by two until a bracket is found,
    then bisect it. amax limits the largest absolute log-velocity update
    to max_step_lnvs. Invalid models form the high end of a bracket.
    At most max_backtracks+1 trial evaluations are allowed per direction.
    At the step cap, also accept an Armijo-decreasing trial whose slope is
    still negative, even if strong curvature is not satisfied. Such a capped
    step is logged explicitly; it is not a strong Wolfe step. Overshooting
    or invalid trials still shrink the bracket. Return None if the trial
    budget or numerical precision prevents finding an acceptable step.

    :param problem: Inversion data and fixed-layer mask
    :type problem: InversionProblem
    :param controls: Optimizer settings including wolfe_c1 and wolfe_c2
    :type controls: InvPara
    :param forward: Evaluator recomputing the response and kernel per trial
    :type forward: VsappForward
    :param state: Current accepted model, misfit and raw ln(Vs) gradient
    :type state: ForwardState
    :param direction: Descent direction in ln(Vs); fixed entries are ignored
    :type direction: numpy.ndarray
    :return: Accepted matching model/kernel state, or None on failure
    :rtype: ForwardState or None
    """
    direction = direction.copy()
    direction[problem.fixed] = 0.
    scale = float(np.max(np.abs(direction)))
    slope0 = float(np.sum(state.gradient * direction))
    if not np.isfinite(scale) or scale == 0 or not np.isfinite(slope0) or slope0 >= 0:
        return None
    logger = getattr(forward, 'logger', None)
    log_speed = np.log(state.model.vs)
    cap = controls.max_step_lnvs / scale
    alpha = min(1., cap)
    low, high = 0., None
    low_misfit = state.misfit
    for _ in range(controls.max_backtracks + 1):
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            speed = np.exp(log_speed + alpha * direction)
        speed[direction == 0] = state.model.vs[direction == 0]
        speed[problem.fixed] = problem.initial_vs[problem.fixed]
        trial = None
        if np.isfinite(speed).all() and np.all(speed > 0):
            try:
                # Recompute synthetic Vsapp and its kernel at this trial model.
                trial = forward.evaluate(speed)
            except ValueError as error:
                if logger is not None:
                    logger.debug('Rejected Wolfe trial alpha=%g: %s', alpha, error)
        slope = np.nan if trial is None else float(np.sum(trial.gradient * direction))
        valid = trial is not None and np.isfinite(trial.misfit) and np.isfinite(slope)
        armijo = valid and trial.misfit <= state.misfit + controls.wolfe_c1 * alpha * slope0
        curvature = valid and abs(slope) <= -controls.wolfe_c2 * slope0
        if logger is not None:
            logger.debug('Wolfe trial alpha=%g, max_dlnVs=%g, misfit=%g, '
                         'slope=%g, Armijo=%s, curvature=%s',
                         alpha, alpha * scale, np.nan if trial is None else trial.misfit,
                         slope, armijo, curvature)
        if not armijo or (low > 0 and trial.misfit >= low_misfit):
            high = alpha
        else:
            if trial.misfit < state.misfit:
                if curvature:
                    return trial
                # The decreasing function still wants a larger step: use the cap.
                if alpha == cap and slope < 0:
                    if logger is not None:
                        logger.info('Accepted capped L-BFGS step: max_dlnVs=%g; '
                                    'Armijo satisfied, strong curvature not met',
                                    controls.max_step_lnvs)
                    return trial
            if slope >= 0:
                high = alpha
            else:
                low, low_misfit = alpha, trial.misfit
        # The low endpoint has negative slope; the high endpoint either
        # fails sufficient decrease, is invalid, or has nonnegative slope.
        candidate = min(2 * alpha, cap) if high is None else .5 * (low + high)
        if candidate == alpha or candidate == low or candidate == high:
            if logger is not None:
                logger.debug('Wolfe search stopped: no representable next step')
            return None
        alpha = candidate
    if logger is not None:
        logger.debug('Wolfe search exhausted %d trial evaluations', controls.max_backtracks + 1)
    return None
