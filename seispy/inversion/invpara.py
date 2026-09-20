# SPDX-License-Identifier: GPL-3.0-only
"""Vsapp inversion parameters and SeisPy-style INI configuration loading."""

import configparser
from copy import deepcopy
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np

from seispy.vsapp import period_grid, positive, vector

from .preprocess import _kernel_options

_SECTIONS = {
    'data': ('periods_s', 'sigma'),
    'forward': ('rayp', 'f0', 'method', 'dt', 'npts', 'shift', 'pre_filt',
                'itmax', 'minderr', 'wlevel'),
    'inversion': ('optimizer', 'smooth_sigma_km', 'zero_halfspace', 'max_iterations',
                  'max_step_lnvs', 'max_backtracks', 'gradient_tol', 'misfit_tol',
                  'misfit_window', 'lbfgs_memory', 'wolfe_c1', 'wolfe_c2'),
}
_INTEGERS = {'npts', 'itmax', 'max_iterations', 'max_backtracks', 'misfit_window', 'lbfgs_memory'}


@dataclass
class InvPara:
    """All numerical settings for Vsapp inversion, independent of the input data.

    Attributes follow the existing inversion keywords. ``periods_s`` and
    ``sigma`` may be arrays; rayp is in s/km, velocities in km/s, time in s,
    and smooth_sigma_km in km. f0 must be supplied before inversion. rayp
    is required for a direct curve and inferred from retained station events.
    Defaults are shared by keyword, parameter-object and config-file calls.
    max_step_lnvs is the initial SD step or the hard L-BFGS update cap.
    At the cap, L-BFGS accepts an Armijo-decreasing step with negative slope
    even when the strong Wolfe curvature condition is not yet satisfied.
    zero_halfspace masks only the last kernel column; the half-space velocity
    remains free to change through depth smoothing of the search direction.

    The model and observed data remain runtime inputs. Use read_para to load
    [data], [forward] and [inversion] sections, or edit attributes directly.
    Inversion validates and copies the parameters before use.
    """

    periods_s: np.ndarray | None = None
    sigma: float | np.ndarray | None = None
    rayp: float | None = None
    f0: float | None = None
    method: str = 'iter'
    dt: float = .05
    npts: int = 1024
    shift: float = 10.
    pre_filt: tuple | None = None
    itmax: int = 400
    minderr: float = .001
    wlevel: float = .05
    optimizer: str = 'gd'
    smooth_sigma_km: float = .5
    zero_halfspace: bool = True
    max_iterations: int = 100
    max_step_lnvs: float = .05
    max_backtracks: int = 20
    gradient_tol: float = 1e-6
    misfit_tol: float = 1e-8
    misfit_window: int = 5
    lbfgs_memory: int = 10
    wolfe_c1: float = 1e-4
    wolfe_c2: float = .9

    def __post_init__(self):
        self.validate()

    def __str__(self):
        return '\n'.join(f'{item.name}: {getattr(self, item.name)}' for item in fields(self))

    def validate(self):
        """Validate the current settings; missing rayp/f0 are allowed until execution."""
        for name in _INTEGERS:
            value = getattr(self, name)
            minimum = 0 if name in ('max_iterations', 'max_backtracks') else 1
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
                raise TypeError(f'{name} must be an integer')
            if value < minimum:
                raise ValueError(f'{name} must be at least {minimum}')
        if self.npts < 2:
            raise ValueError('npts must be at least 2')
        for name in ('dt', 'max_step_lnvs'):
            setattr(self, name, positive(getattr(self, name), name))
        for name in ('rayp', 'f0'):
            value = getattr(self, name)
            if value is not None:
                setattr(self, name, positive(value, name))
        for name in ('shift', 'minderr', 'wlevel', 'smooth_sigma_km', 'gradient_tol', 'misfit_tol'):
            value = getattr(self, name)
            if not np.isscalar(value) or not np.isfinite(value) or value < 0:
                raise ValueError(f'{name} must be finite and nonnegative')
        if self.optimizer not in ('gd', 'lbfgs'):
            raise ValueError("optimizer must be 'gd' or 'lbfgs'")
        if self.method not in ('iter', 'water'):
            raise ValueError("method must be 'iter' or 'water'")
        if not isinstance(self.zero_halfspace, (bool, np.bool_)):
            raise TypeError('zero_halfspace must be a bool')
        if (not np.isscalar(self.wolfe_c1) or not np.isscalar(self.wolfe_c2)
                or not 0 < self.wolfe_c1 < self.wolfe_c2 < 1):
            raise ValueError('Wolfe constants must satisfy 0 < wolfe_c1 < wolfe_c2 < 1')
        if self.periods_s is not None:
            self.periods_s = period_grid(self.periods_s)
        if self.sigma is not None:
            values = np.asarray(self.sigma, dtype=float)
            if values.ndim > 1 or values.size == 0 or not np.isfinite(values).all():
                raise ValueError('sigma must be a finite scalar or vector')
            if np.any(values <= 0):
                raise ValueError('sigma must be positive')
        if self.pre_filt is not None:
            values = vector(self.pre_filt, 'pre_filt')
            if values.size != 2 or not 0 < values[0] < values[1] < .5 / self.dt:
                raise ValueError('pre_filt needs two increasing frequencies below Nyquist')

    @property
    def kernel_kwargs(self):
        """Return fresh forward options for the selected RF deconvolution method."""
        options = dict(dt=self.dt, npts=self.npts, shift=self.shift, pre_filt=self.pre_filt)
        if self.method == 'iter':
            options.update(itmax=self.itmax, minderr=self.minderr)
        else:
            options['wlevel'] = self.wlevel
        return options

    @classmethod
    def read_para(cls, cfg_file):
        """Read a parameter file using SeisPy's sectioned INI convention.

        Numeric arrays accept commas or whitespace. Blank values retain
        defaults. Unknown sections and keys raise rather than silently
        ignoring misspelled parameters. No Python expressions are evaluated.

        :param cfg_file: INI configuration file
        :type cfg_file: str or pathlib.Path
        :return: Validated inversion parameters
        :rtype: InvPara
        :raises FileNotFoundError: If the file does not exist
        :raises ValueError: If a setting, section or parameter name is invalid
        """
        parser = configparser.ConfigParser(interpolation=None, inline_comment_prefixes=('#', ';'))
        try:
            with Path(cfg_file).expanduser().open(encoding='utf-8') as stream:
                parser.read_file(stream)
        except configparser.Error as error:
            raise ValueError(f'Invalid inversion configuration: {error}') from error
        if parser.defaults():
            raise ValueError('Use explicit [data], [forward] and [inversion] sections')
        parameters = cls()
        for section in parser.sections():
            if section not in _SECTIONS:
                raise ValueError(f'Unknown inversion section: {section}')
            for name, value in parser.items(section):
                if name not in _SECTIONS[section]:
                    raise ValueError(f'Unknown parameter [{section}] {name}')
                if not value.strip():
                    continue
                if name in ('periods_s', 'sigma', 'pre_filt'):
                    values = [float(item) for item in value.replace(',', ' ').split()]
                    value = values[0] if name == 'sigma' and len(values) == 1 else np.array(values)
                elif name in _INTEGERS:
                    value = parser.getint(section, name)
                elif name == 'zero_halfspace':
                    value = parser.getboolean(section, name)
                elif name not in ('method', 'optimizer'):
                    value = parser.getfloat(section, name)
                setattr(parameters, name, value)
        parameters.validate()
        return parameters


def invpara(cfg_file):
    """Read an InvPara object, following the hkpara/ccppara function convention.

    :param cfg_file: INI configuration file
    :type cfg_file: str or pathlib.Path
    :return: Inversion settings
    :rtype: InvPara
    """
    return InvPara.read_para(cfg_file)


def resolve_parameters(para=None, cfg_file=None, **overrides):
    """Own a validated copy; explicit keyword settings override object/file values."""
    if para is not None and cfg_file is not None:
        raise ValueError('Specify para or cfg_file, not both')
    if cfg_file is not None:
        parameters = InvPara.read_para(cfg_file)
    elif para is None:
        parameters = InvPara()
    elif isinstance(para, InvPara):
        parameters = deepcopy(para)
    else:
        raise TypeError('para must be an InvPara')
    kernel_options = overrides.pop('kernel_kwargs', None)
    names = {item.name for item in fields(InvPara)}
    unknown = overrides.keys() - names
    if unknown:
        raise TypeError('Unknown inversion parameters: ' + ', '.join(sorted(unknown)))
    for name, value in overrides.items():
        setattr(parameters, name, deepcopy(value))
    # Preserve the existing kernel_kwargs interface; flat keywords have precedence.
    for name, value in _kernel_options(kernel_options, parameters.method).items():
        if name not in overrides:
            setattr(parameters, name, value)
    parameters.validate()
    return parameters
