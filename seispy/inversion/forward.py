# SPDX-License-Identifier: GPL-3.0-only
"""Recompute coupled material properties, synthetic data and kernel per model."""

from copy import deepcopy
from dataclasses import dataclass

import numpy as np

from seispy.core.depmodel import DepModel
from seispy.vsapp_kernel import VsappKernelResult

from .gradient import vsapp_gradient
from .material import brocher_properties


@dataclass(frozen=True)
class ForwardState:
    """One evaluated model and its matching response, misfit and raw ln(Vs) gradient."""

    model: DepModel
    kernel: VsappKernelResult
    misfit: float
    gradient: np.ndarray
    rms_km_s: float


class VsappForward:
    """Evaluate trial models using fixed data settings and Brocher coupling."""

    def __init__(self, problem, log=None):
        self.problem = problem
        self.logger = log
        self.n_evaluations = 0

    def evaluate(self, speed):
        """Build a new model and kernel, including new chain-rule slopes."""
        self.n_evaluations += 1
        if self.logger is not None:
            self.logger.debug('Forward/kernel evaluation %d', self.n_evaluations)
        problem = self.problem
        # Convert the trial Vs to a complete elastic model.
        model = deepcopy(problem.template)
        model.vs = speed.copy()
        model.vp, model.rho, dvp, drho = brocher_properties(speed)
        # This call computes synthetic RF/Vsapp AND the analytic kernel together.
        kernel = model.vsapp_kernel(
            problem.rayp, problem.periods, problem.f0, method=problem.method,
            zero_halfspace=problem.zero_halfspace, vp_vs_derivative=dvp,
            rho_vs_derivative=drho, **problem.options,
        )
        # Form the data misfit and raw gradient with respect to ln(Vs).
        misfit, gradient = vsapp_gradient(
            kernel, problem.observed, vs_km_s=speed, sigma=problem.uncertainty,
        )
        gradient[problem.fixed] = 0.
        model.model_array = np.column_stack((model.depths_elev, model.vp, model.vs, model.rho))
        model.isrho = True
        rms = float(np.sqrt(np.mean((kernel.vsapp - problem.observed)**2)))
        if self.logger is not None:
            self.logger.debug('Forward misfit=%g, rms=%g km/s', misfit, rms)
        return ForwardState(model, kernel, misfit, gradient, rms)
