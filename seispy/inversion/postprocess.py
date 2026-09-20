# SPDX-License-Identifier: GPL-3.0-only
"""Accepted-model snapshots and public Vsapp inversion results."""

from dataclasses import dataclass, replace

import numpy as np

from seispy.core.depmodel import DepModel
from seispy.vsapp_kernel import VsappKernelResult

from .gradient import smooth_gradient


@dataclass(frozen=True)
class VsappIteration:
    """One accepted model, including iteration zero before any update.

    ``misfit`` is weighted mean squared residual divided by two. ``rms_km_s``
    is unweighted. ``max_update_km_s`` is the largest accepted layer change;
    ``max_update_lnvs`` is the largest absolute change in ln(Vs);
    ``gradient_norm`` is the infinity norm of the raw ln(Vs) gradient.
    ``branch_warning`` records the forward kernel diagnostic. ``vs_km_s``
    is an independent, read-only model snapshot. ``gradient_angle_deg`` and
    ``step_length_lnvs`` record the angle and persistent SD step used to
    reach this model (angle None for the first update). Both are None for
    iteration zero and L-BFGS. The angle compares raw gradients before the update.
    """

    iteration: int
    misfit: float
    rms_km_s: float
    max_update_km_s: float
    max_update_lnvs: float
    gradient_norm: float
    branch_warning: bool
    vs_km_s: np.ndarray
    gradient_angle_deg: float | None = None
    step_length_lnvs: float | None = None



@dataclass(frozen=True)
class VsappInversionResult:
    """Final accepted model, matching kernel, gradients and iteration history.

    ``model`` and ``initial_model`` are independent DepModel copies with
    Brocher-derived Vp and density. ``kernel`` was computed at ``model``.
    ``kernel.jacobian`` retains its derivative with respect to Vs.
    ``gradient`` is the local data derivative with respect to ln(Vs),
    including the chain-rule factor Vs;
    ``smoothed_gradient`` is a diagnostic Gaussian-smoothed data gradient;
    an L-BFGS search direction also depends on its secant history.
    ``converged`` reports numerical stopping, not uniqueness of the model.
    Iteration limits, invalid SD updates and failed L-BFGS line searches
    return converged=False.
    ``event_mask`` identifies retained events for invert_station_vsapp;
    it is None for direct curve inversion. Input objects are never modified.
    """

    model: DepModel
    initial_model: DepModel
    kernel: VsappKernelResult
    periods_s: np.ndarray
    observed_vsapp: np.ndarray
    sigma: np.ndarray
    rayp: float
    f0: float
    gradient: np.ndarray
    smoothed_gradient: np.ndarray
    history: tuple[VsappIteration, ...]
    converged: bool
    message: str
    n_evaluations: int
    event_mask: np.ndarray | None = None
    optimizer: str = 'gd'

    @property
    def predicted_vsapp(self):
        """Return the final synthetic Vsapp curve in km/s."""
        return self.kernel.vsapp

    @property
    def misfit_history(self):
        """Return misfits of accepted models, starting with the initial model."""
        return np.array([item.misfit for item in self.history])

    @property
    def n_iterations(self):
        """Return the number of accepted model updates."""
        return len(self.history) - 1



def record_iteration(state, iteration, update, log_update, gradient_norm,
                     gradient_angle_deg=None, step_length_lnvs=None):
    """Capture one accepted state without retaining mutable model arrays."""
    return VsappIteration(
        iteration, state.misfit, state.rms_km_s, update, log_update, gradient_norm,
        bool(state.kernel.diagnostics.get('branch_warning', False)), _snapshot(state.model.vs),
        gradient_angle_deg, step_length_lnvs,
    )


def build_result(problem, outcome, controls):
    """Assemble the final result from the last accepted optimizer state."""
    state = outcome.state
    smoothed = smooth_gradient(state.gradient, problem.depths_km, controls.smooth_sigma_km,
                               fixed_mask=problem.fixed)
    return VsappInversionResult(
        state.model, outcome.initial_model, state.kernel, _snapshot(problem.periods),
        _snapshot(problem.observed), _snapshot(problem.uncertainty), problem.rayp, problem.f0,
        _snapshot(state.gradient), _snapshot(smoothed), tuple(outcome.history),
        outcome.converged, outcome.message, outcome.n_evaluations, optimizer=controls.optimizer,
    )


def attach_event_mask(result, valid):
    """Attach station QC as an independent, read-only snapshot."""
    return replace(result, event_mask=_snapshot(valid))


def _snapshot(values):
    """Create an independent read-only array for returned data and history."""
    result = np.array(values, copy=True)
    result.setflags(write=False)
    return result

