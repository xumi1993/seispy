# SPDX-License-Identifier: GPL-3.0-only
"""Limited-memory inverse-Hessian directions using raw ln(Vs) gradients."""

from collections import deque

import numpy as np


class LBFGSHistory:
    """Store positive-curvature secant pairs and apply the two-loop recursion."""

    def __init__(self, memory):
        self.pairs = deque(maxlen=memory)

    def clear(self):
        """Discard the inverse-Hessian approximation before a gradient restart."""
        self.pairs.clear()

    def update(self, previous, current, fixed):
        """Add s=delta ln(Vs), y=delta raw gradient after an accepted step.

        Fixed coordinates are excluded. Skip nonpositive or numerically
        unreliable curvature rather than changing y or using smoothed gradients.
        """
        step = np.log(current.model.vs) - np.log(previous.model.vs)
        change = current.gradient - previous.gradient
        step[fixed], change[fixed] = 0., 0.
        curvature = float(np.sum(step * change))
        size = float(np.linalg.norm(step) * np.linalg.norm(change))
        if not np.isfinite(curvature) or not np.isfinite(size) or curvature <= 1e-10 * size:
            return False
        self.pairs.append((step.copy(), change.copy(), 1. / curvature))
        return True

    def direction(self, gradient):
        """Return -H*g via the L-BFGS two-loop recursion in log coordinates."""
        q = gradient.copy()
        coefficients = []
        for step, change, reciprocal in reversed(self.pairs):
            alpha = reciprocal * float(np.sum(step * q))
            coefficients.append(alpha)
            q -= alpha * change
        if self.pairs:
            step, change, _ = self.pairs[-1]
            q *= float(np.sum(step * change) / np.sum(change * change))
        for (step, change, reciprocal), alpha in zip(
                self.pairs, reversed(coefficients), strict=True):
            beta = reciprocal * float(np.sum(change * q))
            q += step * (alpha - beta)
        return -q

