"""Independent SynSeis/deconvolution checks for the analytic layered kernels.

Reference waveforms and perturbed models use the original forward solver and
deconvolution functions, never the analytic kernel's forward implementation.
"""

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from seispy import DepModel, SynSeis, decon
from seispy.vsapp import _compute_vsapp_from_components as compute_vsapp
from seispy.vsapp_kernel import (
    deconit_tangent,
    forward_and_jacobian,
    forward_and_jacobian_iter,
    vsapp_kernel,
)


class _ArgmaxTrace:
    """Observe only deconit's own numpy binding; do not patch global numpy."""

    def __init__(self):
        self.path = []

    def __getattr__(self, name):
        return getattr(np, name)

    def argmax(self, values, *args, **kwargs):
        index = np.argmax(values, *args, **kwargs)
        self.path.append(int(index))
        return index


def _deconit(u, w, dt, **settings):
    trace = _ArgmaxTrace()
    with patch.object(decon, "np", trace):
        rf, rms, last = decon.deconit(u, w, dt, **settings)
    return SimpleNamespace(rf=rf, rms=rms, path=trace.path, last=last)


def _model(nlayers=2):
    if nlayers == 1:
        vs, thickness = np.array([3.5]), np.array([0.0])
    elif nlayers == 2:
        vs, thickness = np.array([1.2, 3.5]), np.array([1.0, 0.0])
    elif nlayers == 8:
        vs = np.array([0.8, 1.15, 1.7, 2.3, 2.8, 3.15, 3.5, 3.8])
        thickness = np.array([0.25, 0.4, 0.65, 1.2, 2.0, 3.0, 5.0, 0.0])
    else:
        raise ValueError("Unsupported test model")
    vp = 2 * vs
    return dict(vp=vp, vs=vs, rho=1.74 * vp**0.25, thickness=thickness)


def _coupling(model, coupled):
    if not coupled:
        return {}
    return dict(vp_vs_derivative=np.full_like(model["vs"], 2.0),
                rho_vs_derivative=0.5 * model["rho"] / model["vp"])


def _perturb(model, column, step, coupled):
    result = {key: value.copy() for key, value in model.items()}
    result["vs"][column] += step
    if coupled:
        # Evaluate the nonlinear relation anew; do not use the supplied slopes.
        result["vp"] = 2 * result["vs"]
        result["rho"] = 1.74 * result["vp"]**0.25
    return result


def _reference(model, method, settings):
    synthetic = SynSeis(SimpleNamespace(**model), settings["p"], settings["dt"],
                        npts=settings["npts"], ipha=1)
    synthetic.run_fwd()
    if settings.get("pre_filt") is not None:
        synthetic.filter(*settings["pre_filt"], order=2, zerophase=True)
    r, z = synthetic.rstream[0].data, synthetic.zstream[0].data
    options = dict(tshift=settings["shift"], f0=settings["f0"])
    if method == "iter":
        options.update(itmax=settings.get("itmax", 400), minderr=settings.get("minderr", 0.001))
        radial = _deconit(r, z, settings["dt"], **options)
        vertical = _deconit(z, z, settings["dt"], **options)
        radial_rf, vertical_rf, path = radial.rf, vertical.rf, radial.path
    else:
        options.update(wlevel=settings.get("wlevel", 0.05), normalize=False)
        radial_rf, _ = decon.deconwater(r, z, settings["dt"], **options)
        vertical_rf, _ = decon.deconwater(z, z, settings["dt"], **options)
        path = None
    times = np.arange(len(radial_rf)) * settings["dt"] - settings["shift"]
    curve = compute_vsapp(times, radial_rf, vertical_rf, rayp=settings["p"],
                          periods_s=settings["periods"])
    if set(curve.status) != {"ok"}:
        raise AssertionError(f"Invalid reference curve: {curve.status}")
    return SimpleNamespace(times=times, radial=radial_rf, vertical=vertical_rf,
                           vsapp=curve.vs_km_s, path=path)


def _finite_difference(model, method, settings, relative_step, coupled=False):
    columns = {name: [] for name in ("radial", "vertical", "vsapp")}
    paths = []
    for j, speed in enumerate(model["vs"]):
        h = relative_step * speed
        plus = _reference(_perturb(model, j, h, coupled), method, settings)
        minus = _reference(_perturb(model, j, -h, coupled), method, settings)
        for name in columns:
            columns[name].append((getattr(plus, name) - getattr(minus, name)) / (2 * h))
        paths.extend((plus.path, minus.path))
    return {name: np.column_stack(values) for name, values in columns.items()}, paths


class LayeredKernelTests(unittest.TestCase):
    def setUp(self):
        # Non-power-of-two length, off-sample P time, and off-grid window edges.
        self.settings = dict(p=0.06, dt=0.05, npts=500, shift=6.013, f0=2.0,
                             periods=np.array([0.1, 0.123, 0.375, 0.975, 1.7, 3.5]))

    def assert_reference_forward(self, actual, reference):
        for name in ("times", "radial", "vertical", "vsapp"):
            np.testing.assert_allclose(getattr(actual, name), getattr(reference, name),
                                       rtol=2e-9, atol=2e-10, err_msg=name)
        if reference.path is not None:
            self.assertEqual(actual.diagnostics["radial_path"], reference.path)

    def assert_reference_derivative(self, actual, finite):
        for output, derivative in (("radial", "dradial_dvs"),
                                   ("vertical", "dvertical_dvs"), ("vsapp", "jacobian")):
            analytic = getattr(actual, derivative)
            np.testing.assert_allclose(analytic, finite[output], rtol=5e-5, atol=5e-7,
                                       err_msg=derivative)
            if np.linalg.norm(analytic) > 1e-8:
                self.assertLess(np.linalg.norm(analytic - finite[output]) /
                                np.linalg.norm(analytic), 1e-5, derivative)

    def test_halfspace_physical_identity(self):
        # Vsapp=Vs and dVsapp/dVs=1 follow independently from the half-space
        # free-surface relation and hold with fixed or constrained Vp/density.
        model = _model(1)
        for method, function in (("water", forward_and_jacobian),
                                 ("iter", forward_and_jacobian_iter)):
            for coupled in (False, True):
                for p in (0.04, 0.08):
                    with self.subTest(method=method, coupled=coupled, p=p):
                        actual = function(**model, **(self.settings | {"p": p}),
                                          **_coupling(model, coupled))
                        np.testing.assert_allclose(actual.vsapp, 3.5, rtol=2e-10, atol=2e-10)
                        np.testing.assert_allclose(actual.jacobian, 1.0, rtol=3e-10, atol=3e-10)

    def test_water_forward_and_gradients_including_prefilter_and_coupling(self):
        for nlayers, coupled, pre_filt in ((2, False, None), (2, True, (0.05, 2.0)),
                                          (8, True, None)):
            model = _model(nlayers)
            settings = self.settings | {"pre_filt": pre_filt}
            with self.subTest(nlayers=nlayers, coupled=coupled, pre_filt=pre_filt):
                actual = forward_and_jacobian(**model, **settings, **_coupling(model, coupled))
                reference = _reference(model, "water", settings)
                self.assertEqual(actual.times.shape, (512,))
                self.assertEqual(actual.jacobian.shape, (len(settings["periods"]), nlayers))
                self.assert_reference_forward(actual, reference)
                finite, _ = _finite_difference(model, "water", settings, 1e-5, coupled)
                self.assert_reference_derivative(actual, finite)

    def test_iterative_forward_and_branch_stable_gradients(self):
        model = _model(2)
        for coupled, pre_filt in ((False, None), (True, None), (False, (0.05, 2.0)),
                                  (True, (0.05, 2.0))):
            settings = self.settings | {"pre_filt": pre_filt}
            with self.subTest(coupled=coupled, pre_filt=pre_filt):
                actual = forward_and_jacobian_iter(**model, **settings,
                                                   **_coupling(model, coupled))
                reference = _reference(model, "iter", settings)
                self.assert_reference_forward(actual, reference)
                # Only accept a step if the original greedy solver actually
                # follows the same radial execution path for every perturbation.
                for h in (1e-5, 1e-6, 1e-7, 1e-8):
                    finite, paths = _finite_difference(model, "iter", settings, h, coupled)
                    if all(path == reference.path for path in paths):
                        break
                else:
                    self.fail("No tested central-difference step preserves the iterative path")
                self.assert_reference_derivative(actual, finite)
                np.testing.assert_allclose(actual.dvertical_dvs, 0.0, atol=3e-11)

    def test_depmodel_entrypoint_against_native_forward(self):
        model = _model(2)
        depmod = DepModel.read_layer_model(np.array([0.0, 1.0]), model["thickness"],
                                          model["vp"], model["vs"], rho=model["rho"])
        for method in ("iter", "water"):
            with self.subTest(method=method):
                settings = self.settings | {"pre_filt": (0.05, 2.0)}
                options = {key: value for key, value in settings.items()
                           if key not in ("p", "periods")}
                actual = vsapp_kernel(depmod, settings["p"], settings["periods"],
                                      method=method, **options)
                self.assert_reference_forward(actual, _reference(model, method, settings))

    def test_inputs_are_not_mutated(self):
        for function in (forward_and_jacobian, forward_and_jacobian_iter):
            model = _model(2)
            before = {key: value.copy() for key, value in model.items()}
            periods = self.settings["periods"].copy()
            function(**model, **self.settings, pre_filt=(0.05, 2.0))
            for key in model:
                np.testing.assert_array_equal(model[key], before[key])
            np.testing.assert_array_equal(self.settings["periods"], periods)

    def test_invalid_model_sampling_and_filter_are_rejected(self):
        model = _model(2)
        for function in (forward_and_jacobian, forward_and_jacobian_iter):
            for changes in ({"vs": [-1.0, 3.5]}, {"vp": [2.4]},
                            {"thickness": [1.0, 1.0]}, {"p": 0.5}, {"npts": 1},
                            {"periods": [0.05]}, {"vp_vs_derivative": [1.0]},
                            {"pre_filt": (2.0, 0.05)}, {"pre_filt": (0.0, 2.0)},
                            {"pre_filt": (0.05, 11.0)}):
                with self.subTest(function=function.__name__, changes=changes):
                    with self.assertRaises(ValueError):
                        function(**(model | self.settings | changes))


class IterativeBranchTests(unittest.TestCase):
    def test_tracing_preserves_native_result_and_numpy(self):
        x = np.arange(256)
        w = np.exp(-((x - 30) / 4)**2)
        u = 0.6 * w + 0.3 * np.roll(w, 13)
        options = dict(tshift=3.013, f0=2.0, itmax=8, minderr=0.0)
        expected = decon.deconit(u, w, 0.05, **options)
        traced = _deconit(u, w, 0.05, **options)
        np.testing.assert_array_equal(traced.rf, expected[0])
        np.testing.assert_array_equal(traced.rms, expected[1])
        self.assertEqual(len(traced.path), expected[2] + 1)
        self.assertIs(decon.np, np)

    def test_low_level_waveform_jacobian_against_native_deconit(self):
        x = np.arange(237)
        w = np.exp(-((x - 30) / 4)**2)
        u = 0.6 * w + 0.3 * np.roll(w, 13) - 0.2 * np.roll(w, 37)
        du = np.column_stack((np.roll(w, 13), np.roll(w, 37)))
        dw = np.column_stack((0.1 * np.roll(w, 2), -0.2 * np.roll(w, 3)))
        options = dict(tshift=3.013, f0=2.0, itmax=8, minderr=0.0)
        actual = deconit_tangent(u, w, du, dw, 0.05, **options)
        baseline = _deconit(u, w, 0.05, **options)
        self.assertEqual(actual.path, baseline.path)
        np.testing.assert_allclose(actual.rf, baseline.rf, rtol=1e-12, atol=1e-13)
        h = 1e-6
        for j in range(du.shape[1]):
            plus = _deconit(u + h * du[:, j], w + h * dw[:, j], 0.05, **options)
            minus = _deconit(u - h * du[:, j], w - h * dw[:, j], 0.05, **options)
            self.assertEqual(plus.path, baseline.path)
            self.assertEqual(minus.path, baseline.path)
            np.testing.assert_allclose(actual.jacobian[:, j], (plus.rf - minus.rf) / (2 * h),
                                       rtol=3e-6, atol=3e-9)

    def test_self_deconvolution_has_zero_tangent(self):
        x = np.arange(256)
        w = np.exp(-((x - 30) / 4)**2) + 0.2 * np.exp(-((x - 60) / 7)**2)
        dw = np.column_stack((np.roll(w, 5), np.roll(w, 19)))
        actual = deconit_tangent(w, w, dw, dw, 0.05, tshift=3.0)
        np.testing.assert_allclose(actual.jacobian, 0.0, atol=2e-12)
        np.testing.assert_allclose(actual.rf, decon.deconit(w, w, 0.05, tshift=3.0)[0],
                                   rtol=2e-12, atol=2e-13)

    def test_equal_peak_switch_has_no_smooth_gradient(self):
        w, u, du = np.zeros(64), np.zeros(64), np.zeros((64, 1))
        w[0], u[4], u[12], du[4, 0] = 1.0, 1.0, 1.0, 1.0
        options = dict(tshift=0.0, f0=1e6, itmax=1, minderr=0.0)
        actual = deconit_tangent(u, w, du, np.zeros_like(du), 0.1, **options)
        norms = []
        for h in (1e-4, 1e-5):
            plus = _deconit(u + h * du[:, 0], w, 0.1, **options)
            minus = _deconit(u - h * du[:, 0], w, 0.1, **options)
            self.assertNotEqual(plus.path, minus.path)
            norms.append(np.linalg.norm((plus.rf - minus.rf) / (2 * h)))
        self.assertGreater(norms[1] / norms[0], 9.9)
        self.assertGreater(norms[1], 1e3 * np.linalg.norm(actual.jacobian))

    def test_stopping_threshold_switch_changes_iteration_count(self):
        w, u = np.zeros(64), np.zeros(64)
        w[0], u[4], u[12] = 1.0, 1.0, 0.8
        # With an almost-flat Gaussian the first fitted spike removes energy 1
        # out of 1+0.8^2, so this threshold bisects the two perturbed fits.
        threshold = 100.0 / (1.0 + 0.8**2)
        options = dict(tshift=0.0, f0=1e6, itmax=2, minderr=threshold)
        plus, minus = u.copy(), u.copy()
        plus[12] += 1e-5
        minus[12] -= 1e-5
        first, second = _deconit(plus, w, 0.1, **options), _deconit(minus, w, 0.1, **options)
        self.assertEqual({len(first.path), len(second.path)}, {1, 2})
        self.assertGreater(np.linalg.norm((first.rf - second.rf) / 2e-5), 1e4)


if __name__ == "__main__":
    unittest.main()
