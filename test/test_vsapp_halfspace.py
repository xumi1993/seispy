"""Half-space masking acts on inversion derivatives, preserving the forward."""

import unittest

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from seispy import DepModel, SynSeis
from seispy.vsapp_kernel import forward_and_jacobian, forward_and_jacobian_iter


class TestZeroHalfspace(unittest.TestCase):
    def setUp(self):
        self.model = dict(vp=np.array([2.4, 5.5, 7.]), vs=np.array([1.2, 2.8, 3.5]),
                          rho=np.array([2.1, 2.5, 2.8]), thickness=np.array([1., 2., 0.]))
        self.settings = dict(p=.06, periods=np.array([.2, .5, 2.]), f0=3.5,
                             dt=.05, npts=512, shift=5.)

    def tearDown(self):
        plt.close('all')

    def assert_masked(self, base, masked):
        for name in ('times', 'radial', 'vertical', 'vsapp', 'periods_s', 'thickness_km'):
            np.testing.assert_array_equal(getattr(masked, name), getattr(base, name))
        for name in ('jacobian', 'dradial_dvs', 'dvertical_dvs'):
            before, after = getattr(base, name), getattr(masked, name)
            self.assertEqual(before.shape, after.shape)
            np.testing.assert_array_equal(after[:, :-1], before[:, :-1])
            np.testing.assert_array_equal(after[:, -1], 0.)
        self.assertFalse(base.diagnostics['zero_halfspace'])
        self.assertTrue(masked.diagnostics['zero_halfspace'])

    def test_both_methods_preserve_forward_and_other_derivatives_with_coupling(self):
        for function in (forward_and_jacobian_iter, forward_and_jacobian):
            for coupled in (False, True):
                with self.subTest(function=function.__name__, coupled=coupled):
                    slopes = (dict(vp_vs_derivative=np.full(3, 1.8),
                                   rho_vs_derivative=np.full(3, .1)) if coupled else {})
                    snapshots = {name: value.copy() for name, value in self.model.items()}
                    base = function(**self.model, **self.settings, **slopes)
                    masked = function(**self.model, **self.settings, **slopes,
                                      zero_halfspace=True)
                    self.assertGreater(np.max(np.abs(base.jacobian[:, -1])), 0.)
                    self.assert_masked(base, masked)
                    for name, expected in snapshots.items():
                        np.testing.assert_array_equal(self.model[name], expected)

    def test_model_and_synthetic_adapters_mask_every_ray_parameter(self):
        model = DepModel.read_layer_model(np.array([0., 1., 3.]), self.model['thickness'],
                                         self.model['vp'], self.model['vs'], self.model['rho'])
        before = model.vs.copy()
        for method in ('iter', 'water'):
            with self.subTest(method=method):
                options = dict(method=method, f0=3.5, dt=.05, npts=512, shift=5.)
                base = model.vsapp_kernel(.06, self.settings['periods'], **options)
                masked = model.vsapp_kernel(.06, self.settings['periods'],
                                            zero_halfspace=True, **options)
                self.assert_masked(base, masked)
                synthetic = SynSeis(model, [.05, .07], .05, npts=512)
                synth_options = dict(method=method, f0=3.5, shift=5.)
                unmasked = synthetic.vsapp_kernel(self.settings['periods'], **synth_options)
                masked_rays = synthetic.vsapp_kernel(self.settings['periods'],
                                                     zero_halfspace=True, **synth_options)
                self.assertEqual(len(masked_rays), 2)
                for original, result in zip(unmasked, masked_rays, strict=True):
                    self.assert_masked(original, result)
        np.testing.assert_array_equal(model.vs, before)

    def test_mask_applies_to_gradient_and_colormap_not_only_display(self):
        masked = forward_and_jacobian_iter(**self.model, **self.settings, zero_halfspace=True)
        residual = np.array([1., -2., .5])
        gradient = np.einsum('ij,i->j', masked.jacobian, residual)
        self.assertEqual(gradient[-1], 0.)
        _, ax = masked.plot()
        displayed = np.asarray(ax.collections[0].get_array()).reshape(masked.jacobian.T.shape)
        np.testing.assert_array_equal(displayed[-1], 0.)
        limit = np.max(np.abs(masked.jacobian[:, :-1]))
        np.testing.assert_allclose(ax.collections[0].get_clim(), [-limit, limit])

    def test_single_halfspace_masks_the_only_column_without_changing_vsapp(self):
        model = dict(vp=[7.], vs=[3.5], rho=[2.8], thickness=[0.])
        for function in (forward_and_jacobian_iter, forward_and_jacobian):
            with self.subTest(function=function.__name__):
                result = function(**model, **self.settings, zero_halfspace=True)
                np.testing.assert_allclose(result.vsapp, 3.5)
                for name in ('jacobian', 'dradial_dvs', 'dvertical_dvs'):
                    np.testing.assert_array_equal(getattr(result, name), 0.)

    def test_option_accepts_boolean_only(self):
        for function in (forward_and_jacobian_iter, forward_and_jacobian):
            for value in ('False', None, 0, 1):
                with self.subTest(function=function.__name__, value=value):
                    with self.assertRaisesRegex(TypeError, 'zero_halfspace must be a bool'):
                        function(**self.model, **self.settings, zero_halfspace=value)


if __name__ == '__main__':
    unittest.main()
