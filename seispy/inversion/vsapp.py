# SPDX-License-Identifier: GPL-3.0-only
"""Public entry points coordinating preprocessing, optimization and results."""

from copy import deepcopy

import numpy as np

from seispy.setuplog import inversion_logger

from .forward import VsappForward
from .gradient import smooth_gradient
from .invpara import resolve_parameters
from .lbfgs import LBFGSHistory
from .linesearch import strong_wolfe
from .optimize import (
    OptimizationOutcome,
    gradient_angle,
    gradient_norm,
    windowed_misfit_change,
)
from .postprocess import attach_event_mask, build_result, record_iteration
from .preprocess import prepare_inversion, prepare_station_data


def invert_vsapp(initial_model, periods_s=None, observed_vsapp=None, *,
                 para=None, cfg_file=None, log=None, **kwargs):
    """Invert one Vsapp curve using gradient descent or L-BFGS in ln(Vs).

    At every trial model, recompute Vp(Vs), density(Vp) and their analytic
    slopes using Brocher (2005), then call DepModel.vsapp_kernel. Thus every
    accepted update has a newly computed kernel; no initial kernel is reused
    after changing the model. Thickness, ray parameter and RF processing
    settings stay fixed. Initial Vp and density are replaced in an owned copy
    by the same Brocher relations used in all later models.

    The optimization variable is ln(Vs), with gradient Vs * d(misfit)/d(Vs).
    For optimizer='gd' (steepest descent), normalize the depth-smoothed
    negative gradient by its largest absolute entry. Start with a log step
    of max_step_lnvs. Compare raw gradients at successive accepted models;
    an angle strictly greater than 120 degrees halves the persistent step.
    Allow 1e-12 degrees of roundoff at the threshold. Otherwise retain it.
    Update Vs_new = Vs * exp(delta_lnVs) and evaluate the new model once.
    There is no Armijo check, backtracking or velocity clipping in SD mode.
    Misfit may increase. Convergence compares absolute changes of window means,
    so a small increase or an oscillation with stable means may satisfy it.
    An invalid model stops SD at the last valid state without another trial.

    With optimizer='lbfgs', use limited-memory BFGS secant pairs in ln(Vs)
    and a strong Wolfe line search (Armijo plus absolute curvature condition).
    Both secant pairs and Wolfe derivatives use raw gradients. Gaussian
    smoothing acts on the proposed L-BFGS direction; retry without smoothing,
    then restart with -g if needed. Bracketing can expand a step up to
    max_step_lnvs, then zoom bisects the bracket. At this hard update cap,
    accept a strictly decreasing Armijo step with negative directional
    derivative even if the strong curvature condition is not met. Log this
    capped-step fallback explicitly. Overshooting or invalid models still
    require backtracking. Secant pairs with unreliable curvature are skipped.

    With zero_halfspace=True, only the last kernel column is zeroed. The
    half-space still participates in depth smoothing and velocity updates,
    so its update can come from neighbouring layers despite its zero raw
    gradient. With iterative deconvolution the kernel is local to an execution
    branch; branch warnings are retained in history and line-search failure
    is reported explicitly.

    Parameters are centralized in InvPara. Supply para=InvPara(...),
    cfg_file='inversion.cfg', or individual keyword overrides. Explicit
    keywords take precedence over the file/object; periods_s supplied as
    an argument overrides its configured value. Input parameter objects
    are copied. Model and observed_vsapp remain runtime data inputs.
    Progress and termination use SetupLog.Invlog or the supplied logger.

    :param para: Inversion parameter object, defaults to None
    :type para: InvPara, optional
    :param cfg_file: INI file, mutually exclusive with para, defaults to None
    :type cfg_file: str or pathlib.Path, optional
    :param log: Existing SeisPy SetupLog or Python logger, defaults to None
    :type log: object, optional
    :param kwargs: InvPara keyword overrides; legacy kernel_kwargs is also accepted
    :type kwargs: dict
    :param initial_model: Initial sampled velocity model; input is not modified
    :type initial_model: seispy.DepModel
    :param periods_s: Increasing cosine-squared smoothing half-widths in seconds
    :type periods_s: array_like
    :param observed_vsapp: Positive observed curve in km/s, one value per period
    :type observed_vsapp: array_like
    :param rayp: Fixed P-wave ray parameter in s/km; use retained-event mean
        slowness when approximating a station mean curve
    :type rayp: float
    :param f0: Gaussian factor used for RF deconvolution
    :type f0: float
    :param sigma: Positive scalar or per-period uncertainty in km/s.
        None uses equal unit weights, defaults to None
    :type sigma: float or array_like, optional
    :param method: RF deconvolution, 'iter' or 'water', defaults to 'iter'
    :type method: str, optional
    :param smooth_sigma_km: Gaussian standard deviation along depth in km.
        Zero disables smoothing, defaults to 0.5
    :type smooth_sigma_km: float, optional
    :param zero_halfspace: Zero the half-space kernel column without fixing its
        velocity; depth smoothing can still update it, defaults to True
    :type zero_halfspace: bool, optional
    :param max_iterations: Maximum number of accepted updates, defaults to 100
    :type max_iterations: int, optional
    :param max_step_lnvs: Initial persistent SD log step, or maximum absolute
        log update in L-BFGS, dimensionless. Default 0.05 permits factors
        exp(-0.05) to exp(0.05), approximately -4.88% to +5.13%
    :type max_step_lnvs: float, optional
    :param max_backtracks: L-BFGS allows max_backtracks+1 trials per direction
        across bracketing and zoom. Unused by SD, defaults to 20
    :type max_backtracks: int, optional
    :param gradient_tol: Raw ln(Vs)-gradient infinity-norm stopping tolerance, defaults to 1e-6
    :type gradient_tol: float, optional
    :param misfit_tol: Stop when abs(A-B)/B < misfit_tol, where A and B are the
        mean weighted misfits in the last and preceding misfit_window accepted
        updates. Requires at least 2*misfit_window updates; iteration zero is
        excluded. Zero disables this criterion, defaults to 1e-8
    :type misfit_tol: float, optional
    :param misfit_window: Positive number of updates in each of two adjacent,
        nonoverlapping misfit windows, defaults to 5. Both zero means give
        zero relative change; B=0 with A>0 never passes this criterion.
        Gradient tolerance and iteration limits are checked independently
    :type misfit_window: int, optional
    :param optimizer: 'gd' for gradient-angle-controlled steepest descent or 'lbfgs' for strong
        Wolfe L-BFGS, defaults to 'gd'; independent of RF deconvolution method
    :type optimizer: str, optional
    :param lbfgs_memory: Maximum retained secant pairs, defaults to 10
    :type lbfgs_memory: int, optional
    :param wolfe_c1: Sufficient-decrease coefficient, defaults to 1e-4
    :type wolfe_c1: float, optional
    :param wolfe_c2: Strong-curvature coefficient, defaults to 0.9.
        Must satisfy 0 < wolfe_c1 < wolfe_c2 < 1
    :type wolfe_c2: float, optional
    :param kernel_kwargs: Fixed dt, npts, shift, pre_filt and method-specific
        itmax/minderr or wlevel options, defaults to None
    :type kernel_kwargs: dict, optional
    :return: Final accepted model, synthetic curve, kernel and convergence history
    :rtype: VsappInversionResult
    :raises ValueError: If inputs, initial model or initial forward calculation are invalid
    :raises TypeError: If model, options or option types are invalid
    """
    # 1. Load/copy all settings, then prepare independent data and model copies.
    logger = inversion_logger(log)
    if periods_s is not None:
        kwargs['periods_s'] = periods_s
    try:
        controls = resolve_parameters(para, cfg_file, **kwargs)
        if controls.f0 is None or controls.rayp is None:
            raise ValueError('f0 and rayp are required for direct Vsapp inversion')
        if controls.periods_s is None or observed_vsapp is None:
            raise ValueError('periods_s and observed_vsapp are required')
        problem = prepare_inversion(
            initial_model, controls.periods_s, observed_vsapp,
            rayp=controls.rayp, f0=controls.f0, sigma=controls.sigma,
            method=controls.method, zero_halfspace=controls.zero_halfspace,
            kernel_kwargs=controls.kernel_kwargs,
        )
    except (ValueError, TypeError, OSError) as error:
        logger.error('Invalid inversion input: %s', error)
        raise
    logger.info('Starting Vsapp inversion: optimizer=%s, method=%s, periods=%d, layers=%d',
                controls.optimizer, controls.method, problem.periods.size, problem.initial_vs.size)
    logger.debug('Inversion parameters:\n%s', controls)

    # 2. Initial forward calculation: Brocher model -> RF/Vsapp + kernel -> misfit/gradient.
    forward = VsappForward(problem, log=logger)
    try:
        current = forward.evaluate(problem.initial_vs)
    except (ValueError, TypeError) as error:
        logger.error('Initial forward calculation failed: %s', error)
        raise
    logger.info('Iteration 0: misfit=%.8g, rms=%.6g km/s, gradient=%.6g',
                current.misfit, current.rms_km_s, gradient_norm(current))
    start_model = deepcopy(current.model)
    history = [record_iteration(current, 0, 0., 0., gradient_norm(current))]
    converged = False
    message = 'Maximum iterations reached'

    # SD carries its step length across iterations. L-BFGS uses strong Wolfe.
    previous_gradient = None
    sd_step_length = controls.max_step_lnvs
    memory = None
    if controls.optimizer == 'lbfgs':
        memory = LBFGSHistory(controls.lbfgs_memory)

    for iteration in range(1, controls.max_iterations + 1):
        # 3. The current gradient belongs to the current model, not the initial model.
        if gradient_norm(current) <= controls.gradient_tol:
            converged = True
            message = 'Log-gradient tolerance reached'
            break

        # 4. Compute a descent direction in ln(Vs), then smooth along depth.
        direction = -current.gradient
        if memory is not None:
            direction = memory.direction(current.gradient)
        direction[problem.fixed] = 0.
        smoothed_direction = smooth_gradient(
            direction, problem.depths_km, controls.smooth_sigma_km,
            fixed_mask=problem.fixed,
        )

        # 5. SD updates once; only an angle >120 degrees halves its persistent step.
        # L-BFGS instead evaluates trial models inside a strong Wolfe line search.
        angle = None
        step_length = None
        if memory is None:
            angle = gradient_angle(previous_gradient, current.gradient)
            if angle is not None and angle > 120. + 1e-12:
                sd_step_length /= 2.
                logger.info('SD gradient angle=%.3f deg >120: step_lnvs reduced to %.6g',
                            angle, sd_step_length)
            step_length = sd_step_length
            scale = float(np.max(np.abs(smoothed_direction)))
            if scale == 0:
                message = 'SD smoothed direction is zero'
                break
            delta_lnvs = sd_step_length * (smoothed_direction / scale)
            with np.errstate(over='ignore', under='ignore', invalid='ignore'):
                speed = np.exp(np.log(current.model.vs) + delta_lnvs)
            speed[smoothed_direction == 0] = current.model.vs[smoothed_direction == 0]
            speed[problem.fixed] = problem.initial_vs[problem.fixed]
            if not np.isfinite(speed).all() or np.any(speed <= 0):
                message = 'SD update produced a nonfinite or nonpositive velocity'
                break
            # Exactly one new forward/kernel evaluation; no Armijo check or backtracking.
            try:
                trial = forward.evaluate(speed)
            except ValueError as error:
                message = f'SD forward calculation failed: {error}'
                break
        else:
            trial = strong_wolfe(problem, controls, forward, current, smoothed_direction)
            if trial is None and not np.array_equal(smoothed_direction, direction):
                logger.debug('Retrying Wolfe search with the unsmoothed L-BFGS direction')
                trial = strong_wolfe(problem, controls, forward, current, direction)
            if trial is None:
                logger.debug('Wolfe search failed; resetting L-BFGS history')
                memory.clear()
                raw_direction = -current.gradient
                if (not np.array_equal(raw_direction, direction)
                        and not np.array_equal(raw_direction, smoothed_direction)):
                    trial = strong_wolfe(problem, controls, forward, current, raw_direction)
            if trial is None:
                message = 'L-BFGS line search failed: no acceptable Wolfe or capped-decrease step'
                break

        # 6. Accept the trial model AND its matching synthetic data/kernel/gradient.
        # Reuse this evaluated state next iteration, with no duplicate forward call.
        if memory is not None:
            memory.update(current, trial, problem.fixed)
        previous_gradient = current.gradient.copy()
        previous = current
        current = trial
        update = float(np.max(np.abs(current.model.vs - previous.model.vs)))
        log_update = float(np.max(np.abs(np.log(current.model.vs) - np.log(previous.model.vs))))
        history.append(record_iteration(current, iteration, update, log_update,
                                         gradient_norm(current), angle, step_length))

        logger.info('Iteration %d: misfit=%.8g, rms=%.6g km/s, gradient=%.6g, dlnVs=%.6g',
                    iteration, current.misfit, current.rms_km_s, gradient_norm(current), log_update)

        # 7. Compare the last N accepted misfits with the preceding N misfits.
        # No window-based stopping until 2N updates; iteration zero is excluded.
        if controls.misfit_tol > 0:
            relative_change = windowed_misfit_change(history, controls.misfit_window)
            if relative_change is not None:
                logger.debug('Window mean misfit change=%g, tolerance=%g',
                             relative_change, controls.misfit_tol)
            if relative_change is not None and relative_change < controls.misfit_tol:
                converged = True
                message = 'Windowed misfit change tolerance reached'
                break

    # Also check the gradient after the final allowed update (or with zero iterations).
    if gradient_norm(current) <= controls.gradient_tol:
        converged = True
        message = 'Log-gradient tolerance reached'
    outcome = OptimizationOutcome(current, start_model, tuple(history), converged,
                                  message, forward.n_evaluations)
    if converged:
        logger.info('Vsapp inversion stopped: %s; iterations=%d, evaluations=%d, misfit=%.8g',
                    message, len(history)-1, forward.n_evaluations, current.misfit)
    else:
        logger.warning('Vsapp inversion stopped: %s; iterations=%d, evaluations=%d, misfit=%.8g',
                       message, len(history)-1, forward.n_evaluations, current.misfit)
    return build_result(problem, outcome, controls)


def invert_station_vsapp(initial_model, measurements, *, para=None, cfg_file=None,
                         log=None, **kwargs):
    """Invert a station mean using the mean ray parameter of retained events.

    Reject a whole event if any period is invalid, as in StationVsAppResult.plot.
    Average the surviving curves and their slownesses, then call invert_vsapp.
    This uses the average-slowness approximation; no moveout is performed.
    No uncertainty is inferred from scatter: pass sigma explicitly to weight
    the mean curve, or leave it unset for equal weights. A common forward f0
    must be supplied to match the observations' processing convention.

    :param initial_model: Initial sampled model
    :type initial_model: seispy.DepModel
    :param measurements: Event-wise station measurements with rayp in s/km
    :type measurements: seispy.vsapp.StationVsAppResult
    :param para: InvPara settings; f0 must be set here, in cfg_file or in kwargs
    :type para: InvPara, optional
    :param cfg_file: INI configuration file, defaults to None
    :type cfg_file: str or pathlib.Path, optional
    :param log: Existing SetupLog or Python logger, defaults to None
    :type log: object, optional
    :param kwargs: Keyword options passed to invert_vsapp, except rayp
    :type kwargs: dict
    :return: Inversion result with an accepted-event mask
    :rtype: VsappInversionResult
    :raises ValueError: If station arrays are invalid or no complete event remains
    :raises TypeError: If measurements has the wrong type or rayp is overridden
    """
    logger = inversion_logger(log)
    try:
        if 'rayp' in kwargs:
            raise TypeError('invert_station_vsapp uses the retained-event mean rayp')
        parameters = resolve_parameters(para, cfg_file, **kwargs)
        periods, observed, rayp, valid = prepare_station_data(measurements)
    except (ValueError, TypeError, OSError) as error:
        logger.error('Invalid station inversion input: %s', error)
        raise
    # Station data determine both the period grid and the mean retained-event ray parameter.
    parameters.periods_s = periods
    parameters.rayp = rayp
    logger.info('Station Vsapp: retained %d/%d events, mean rayp=%.8g s/km',
                np.count_nonzero(valid), valid.size, rayp)
    result = invert_vsapp(initial_model, periods, observed, para=parameters, log=logger)
    return attach_event_mask(result, valid)
