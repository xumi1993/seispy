"""Native SeisPy model and synthetic-waveform integration regressions."""

import unittest
from unittest import mock

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from seispy.core.depmodel import DepModel
from seispy.decon import deconit, deconwater
from seispy.seisfwd import SynSeis
from seispy.utils import vs2vprho
from seispy.vsapp import _compute_vsapp_from_components as compute_vsapp


def layer_model():
    return DepModel.read_layer_model(
        np.array([0., 1.]), np.array([1., 0.]),
        np.array([2.4, 7.]), np.array([1.2, 3.5]),
        rho=np.array([2.1, 2.8]),
    )


class ModelVsAppTests(unittest.TestCase):
    def test_layer_reader_preserves_vp_vs_and_density(self):
        depth = np.arange(8.)
        vp, vs, rho = [4., 6., 8.], [2., 3., 4.], [2., 2.5, 3.]
        model = DepModel.read_layer_model(depth, [3., 2., 0.], vp, vs, rho)
        indices = [0, 0, 0, 1, 1, 2, 2, 2]
        assert_array_equal(model.vp, np.asarray(vp)[indices])
        assert_array_equal(model.vs, np.asarray(vs)[indices])
        assert_array_equal(model.rho, np.asarray(rho)[indices])
        assert_array_equal(model.thickness, [1., 1., 1., 1., 1., 1., 1., 0.])

    def test_layer_reader_density_is_from_supplied_vs(self):
        model = DepModel.read_layer_model(
            np.array([0., 1.]), [1., 0.], [3.8, 6.1], [2., 3.4],
        )
        assert_array_equal(model.vp, [3.8, 6.1])
        assert_array_equal(model.vs, [2., 3.4])
        assert_allclose(model.rho, vs2vprho(np.array([2., 3.4]))[1])

    def test_model_kernel_delegates_explicitly(self):
        model = layer_model()
        periods = np.array([.2, .5])
        with mock.patch('seispy.vsapp_kernel.vsapp_kernel') as wrapped:
            result = model.vsapp_kernel(.06, periods, f0=2., method='water', wlevel=.01)
        self.assertIs(result, wrapped.return_value)
        wrapped.assert_called_once_with(model, .06, periods, f0=2., method='water',
                                        zero_halfspace=False, wlevel=.01)

    def test_synthetic_curves_match_native_deconvolution(self):
        periods = np.array([.2, .5, 1.])
        for method in ('iter', 'water'):
            for pre_filt in (None, (.05, 2.)):
                with self.subTest(method=method, pre_filt=pre_filt):
                    synthetic = SynSeis(layer_model(), [.05, .07], .05, npts=257)
                    options = ({'itmax': 100, 'minderr': .002}
                               if method == 'iter' else {'wlevel': .02})
                    curves = synthetic.compute_vsapp(
                        periods, method=method, shift=5.013, f0=2.,
                        pre_filt=pre_filt, **options,
                    )
                    self.assertIsInstance(curves, tuple)
                    self.assertEqual(len(curves), 2)
                    self.assertFalse(hasattr(synthetic, 'rstream'))
                    synthetic.run_fwd()
                    if pre_filt is not None:
                        synthetic.filter(*pre_filt, order=2, zerophase=True)
                    deconvolve = deconit if method == 'iter' else deconwater
                    for rayp, rtrace, ztrace, curve in zip(
                            synthetic.rayp, synthetic.rstream, synthetic.zstream,
                            curves, strict=True):
                        radial = deconvolve(rtrace.data, ztrace.data, .05,
                                            tshift=5.013, f0=2., **options)[0]
                        vertical = deconvolve(ztrace.data, ztrace.data, .05,
                                              tshift=5.013, f0=2., **options)[0]
                        expected = compute_vsapp(
                            np.arange(len(radial)) * .05 - 5.013, radial, vertical,
                            rayp=rayp, periods_s=periods,
                        )
                        assert_array_equal(curve.vs_km_s, expected.vs_km_s)
                        assert_array_equal(curve.radial_integral, expected.radial_integral)
                        assert_array_equal(curve.vertical_integral, expected.vertical_integral)
                        self.assertEqual(curve.rayp, rayp)

    def test_synthetic_operations_leave_existing_streams_and_model_unchanged(self):
        model = layer_model()
        original_model = [getattr(model, name).copy() for name in ('vp', 'vs', 'rho', 'thickness')]
        synthetic = SynSeis(model, .06, .05, npts=512)
        synthetic.run_fwd()
        synthetic.rstream[0].data[:] = 4.
        synthetic.zstream[0].data[:] = 7.
        rstream, zstream = synthetic.rstream, synthetic.zstream
        synthetic.compute_vsapp([.2, .5], shift=5., pre_filt=(.05, 2.))
        synthetic.vsapp_kernel([.2, .5], shift=5., pre_filt=(.05, 2.))
        self.assertIs(synthetic.rstream, rstream)
        self.assertIs(synthetic.zstream, zstream)
        assert_array_equal(synthetic.rstream[0].data, np.full(512, 4.))
        assert_array_equal(synthetic.zstream[0].data, np.full(512, 7.))
        for name, expected in zip(('vp', 'vs', 'rho', 'thickness'), original_model,
                                  strict=True):
            assert_array_equal(getattr(model, name), expected)

    def test_synthetic_kernel_matches_native_vsapp(self):
        periods = np.array([.2, .5, 1.])
        synthetic = SynSeis(layer_model(), [.05, .07], .05, npts=512)
        for method in ('iter', 'water'):
            for pre_filt in (None, (.05, 2.)):
                with self.subTest(method=method, pre_filt=pre_filt):
                    options = dict(method=method, shift=5., f0=2., pre_filt=pre_filt)
                    measured = synthetic.compute_vsapp(periods, **options)
                    kernels = synthetic.vsapp_kernel(periods, **options)
                    self.assertIsInstance(kernels, tuple)
                    self.assertEqual(len(kernels), 2)
                    for curve, kernel in zip(measured, kernels, strict=True):
                        assert_allclose(kernel.vsapp, curve.vs_km_s, rtol=1e-8, atol=1e-9)
                        self.assertEqual(kernel.jacobian.shape, (3, 2))
                        self.assertTrue(np.isfinite(kernel.jacobian).all())

    def test_kernel_preserves_sampled_model_parameter_count(self):
        model = DepModel.read_layer_model(
            np.arange(4.), [2., 0.], [4., 7.], [2., 3.5], rho=[2.3, 2.8],
        )
        kernel = model.vsapp_kernel(.06, [.2, .5], f0=2., dt=.05, npts=512, shift=5.)
        self.assertEqual(kernel.jacobian.shape, (2, len(model.vs)))
        self.assertEqual(kernel.jacobian.shape[1], 4)

    def test_synthetic_rejects_non_p_incidence(self):
        synthetic = SynSeis(layer_model(), .06, .05, ipha=-1)
        for method in (synthetic.compute_vsapp, synthetic.vsapp_kernel):
            with self.assertRaisesRegex(ValueError, 'ipha=1'):
                method([.2])

    def test_synthetic_rejects_unknown_method_and_parameter_overrides(self):
        synthetic = SynSeis(layer_model(), .06, .05)
        for method in (synthetic.compute_vsapp, synthetic.vsapp_kernel):
            with self.assertRaises(ValueError):
                method([.2], method='unknown')
            for name in ('dt', 'npts', 'p', 'phase', 'normalize', 'ipha', 'depmod'):
                with self.subTest(method=method.__name__, option=name):
                    with self.assertRaises(TypeError):
                        method([.2], **{name: 1})

    def test_synthetic_rejects_invalid_bandpass(self):
        synthetic = SynSeis(layer_model(), .06, .05)
        for pre_filt in ((0., 2.), (2., 1.), (1., 10.), (1., np.nan), (1.,)):
            with self.subTest(pre_filt=pre_filt):
                with self.assertRaises(ValueError):
                    synthetic.compute_vsapp([.2], pre_filt=pre_filt)

    def test_synthetic_rejects_model_length_mismatch_before_native_forward(self):
        model = layer_model()
        model.thickness = np.array([1., 1., 0.])
        synthetic = SynSeis(model, .06, .05)
        with mock.patch.object(SynSeis, 'run_fwd') as forward:
            for method in (synthetic.compute_vsapp, synthetic.vsapp_kernel):
                with self.assertRaisesRegex(ValueError, 'equal length'):
                    method([.2])
        forward.assert_not_called()

    def test_synthetic_rejects_invalid_sampling_and_ray_parameters(self):
        cases = (
            {'dt': 0.}, {'dt': np.nan}, {'npts': 2}, {'npts': 512.5},
            {'rayp': []}, {'rayp': [-.06]}, {'rayp': [np.nan]},
            {'rayp': [[.06]]}, {'rayp': [.5]},
        )
        for kwargs in cases:
            with self.subTest(kwargs=kwargs):
                parameters = dict(rayp=.06, dt=.05, npts=512)
                parameters.update(kwargs)
                synthetic = SynSeis(layer_model(), **parameters)
                with mock.patch.object(SynSeis, 'run_fwd') as forward:
                    for method in (synthetic.compute_vsapp, synthetic.vsapp_kernel):
                        with self.assertRaises(ValueError):
                            method([.2])
                forward.assert_not_called()


if __name__ == '__main__':
    unittest.main()
