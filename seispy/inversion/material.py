# SPDX-License-Identifier: GPL-3.0-only
"""Brocher (2005) material relations and their analytic chain-rule slopes."""

import numpy as np

from seispy.utils import vs2vprho
from seispy.vsapp import vector


def brocher_properties(vs):
    """Return Vp, density and their total derivatives with respect to Vs.

    Values reuse :func:`seispy.utils.vs2vprho`: the Vp(Vs) quartic and the
    density(Vp) quintic in Brocher (2005), doi:10.1785/0120050077. The
    density derivative includes dVp/dVs via the chain rule. Velocities are
    in km/s and densities in g/cm^3. No clipping of the relations is applied.

    :param vs: Positive layer S velocities in km/s
    :type vs: array_like
    :return: Vp, density, dVp/dVs and dDensity/dVs arrays
    :rtype: (numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray)
    :raises ValueError: If velocities or derived physical properties are invalid
    """
    speed = vector(vs, 'vs')
    if np.any(speed <= 0):
        raise ValueError('vs must be positive')
    vp, rho = vs2vprho(speed)
    dvp = 2.0947 - 2 * .8206 * speed + 3 * .2683 * speed**2 - 4 * .0251 * speed**3
    drho_dvp = (1.6612 - 2 * .4721 * vp + 3 * .0671 * vp**2
                - 4 * .0043 * vp**3 + 5 * .000106 * vp**4)
    drho = drho_dvp * dvp
    if (not all(np.isfinite(a).all() for a in (vp, rho, dvp, drho))
            or np.any(vp <= speed) or np.any(rho <= 0)):
        raise ValueError('Brocher relations require finite Vp > Vs and positive density')
    return vp, rho, dvp, drho
