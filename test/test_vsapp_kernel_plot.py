"""Kernel plots preserve physical coordinates and signed layer derivatives."""

import unittest
from dataclasses import replace
from unittest.mock import patch

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from seispy import DepModel, SynSeis
from seispy.vsapp_kernel import (
    VsappKernelResult,
    forward_and_jacobian,
    forward_and_jacobian_iter,
)


class TestVsappKernelPlot(unittest.TestCase):
    def setUp(self):
        self.result = VsappKernelResult(
            times=np.array([0.]), radial=np.array([1.]), vertical=np.array([1.]),
            vsapp=np.array([2., 2.5, 3.]), dradial_dvs=np.zeros((1, 3)),
            dvertical_dvs=np.zeros((1, 3)),
            jacobian=np.array([[1., -.5, .2], [.7, .1, -.1], [.2, .8, .3]]),
            diagnostics={}, periods_s=np.array([.2, .6, 1.4]),
            thickness_km=np.array([1., 3., 0.]),
        )

    def tearDown(self):
        plt.close('all')

    def test_nonuniform_periods_and_layers_display_the_correct_derivatives(self):
        fig, ax = self.result.plot()
        mesh = ax.collections[0]
        np.testing.assert_array_equal(np.asarray(mesh.get_array()).reshape(3, 3),
                                      self.result.jacobian.T)
        coordinates = mesh.get_coordinates()
        np.testing.assert_allclose(coordinates[0, :, 0], [0., .4, 1., 1.8], atol=1e-14)
        np.testing.assert_allclose(coordinates[:, 0, 1], [0., 1., 4., 7.])
        self.assertEqual(mesh.get_clim(), (-1., 1.))
        np.testing.assert_allclose(ax.get_ylim(), [7., 0.])
        self.assertEqual(len(fig.axes), 2)
        self.assertIn('Half-space', [text.get_text() for text in ax.texts])

    def test_existing_axes_shared_color_limits_and_depth_clipping(self):
        fig, ax = plt.subplots()
        returned_fig, returned_ax = self.result.plot(
            ax=ax, max_depth=3., vmax=.5, colorbar=False, cmap='PiYG', title='Test',
        )
        self.assertIs(returned_fig, fig)
        self.assertIs(returned_ax, ax)
        self.assertEqual(len(fig.axes), 1)
        self.assertEqual(ax.collections[0].get_clim(), (-.5, .5))
        self.assertEqual(ax.collections[0].get_cmap().name, 'PiYG')
        self.assertEqual(ax.get_title(), 'Test')
        np.testing.assert_allclose(ax.get_ylim(), [3., 0.])
        self.assertEqual(len(ax.texts), 0)

    def test_single_period_homogeneous_halfspace_and_zero_kernel(self):
        result = replace(self.result, periods_s=np.array([.4]),
                         thickness_km=np.array([0.]), jacobian=np.zeros((1, 1)))
        for bottom in (None, 5.):
            with self.subTest(max_depth=bottom):
                _, ax = result.plot(max_depth=bottom)
                mesh = ax.collections[0]
                self.assertEqual(mesh.get_clim(), (-1., 1.))
                np.testing.assert_allclose(mesh.get_coordinates()[0, :, 0], [.2, .6])
                np.testing.assert_allclose(ax.get_ylim(), [bottom or 1., 0.])
                self.assertEqual(ax.texts[0].get_text(), 'Half-space')

    def test_plot_does_not_modify_results_or_display_unless_requested(self):
        snapshots = {name: value.copy() for name, value in vars(self.result).items()
                     if isinstance(value, np.ndarray)}
        with patch('matplotlib.pyplot.show') as show:
            self.result.plot()
            show.assert_not_called()
            self.result.plot(show=True)
            show.assert_called_once_with()
        for name, expected in snapshots.items():
            np.testing.assert_array_equal(getattr(self.result, name), expected)

    def test_invalid_coordinates_and_limits_fail_before_creating_a_figure(self):
        cases = [
            (replace(self.result, jacobian=np.zeros((2, 3))), {}),
            (replace(self.result, jacobian=np.full((3, 3), np.nan)), {}),
            (replace(self.result, periods_s=np.array([.2, .1, 1.])), {}),
            (replace(self.result, thickness_km=np.array([1., 0., 0.])), {}),
            (replace(self.result, thickness_km=np.array([1., 3., 2.])), {}),
        ]
        cases += [(self.result, {name: value}) for name in ('vmax', 'max_depth')
                  for value in (0., -1., np.nan, np.inf)]
        for result, options in cases:
            with self.subTest(options=options):
                count = len(plt.get_fignums())
                with self.assertRaises(ValueError):
                    result.plot(**options)
                self.assertEqual(len(plt.get_fignums()), count)

    def test_all_forward_entrypoints_return_named_results_with_own_coordinates(self):
        for method, function in (('iter', forward_and_jacobian_iter),
                                 ('water', forward_and_jacobian)):
            with self.subTest(method=method):
                periods = np.array([.2, .5, 1.])
                thickness = np.array([1., 0.])
                vp, vs, rho = np.array([2.4, 7.]), np.array([1.2, 3.5]), np.array([2.1, 2.8])
                model = DepModel.read_layer_model(np.array([0., 1.]), thickness, vp, vs, rho)
                direct = function(vp, vs, rho, thickness, periods=periods,
                                  dt=.05, npts=512, shift=5.)
                adapted = model.vsapp_kernel(.06, periods, f0=2., method=method,
                                            dt=.05, npts=512, shift=5.)
                synthetic = SynSeis(model, .06, .05, npts=512).vsapp_kernel(
                    periods, method=method, shift=5.,
                )[0]
                for result in (direct, adapted, synthetic):
                    self.assertIsInstance(result, VsappKernelResult)
                    np.testing.assert_allclose(result.jacobian, direct.jacobian)
                    np.testing.assert_array_equal(result.periods_s, periods)
                    np.testing.assert_array_equal(result.thickness_km, thickness)
                    self.assertFalse(np.shares_memory(result.periods_s, periods))
                    self.assertFalse(np.shares_memory(result.thickness_km, model.thickness))
                    result.plot()
                periods[:] = 5.
                thickness[:] = 2.
                model.thickness[:] = 4.
                np.testing.assert_array_equal(direct.periods_s, [.2, .5, 1.])
                np.testing.assert_array_equal(direct.thickness_km, [1., 0.])
                np.testing.assert_array_equal(adapted.thickness_km, [1., 0.])


if __name__ == '__main__':
    unittest.main()
