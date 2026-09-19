"""Station plots summarize valid event measurements without changing them."""

import unittest
from dataclasses import replace
from unittest.mock import patch

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from seispy.vsapp import StationVsAppResult, VsAppError


class TestStationVsAppPlot(unittest.TestCase):
    def setUp(self):
        self.result = StationVsAppResult(
            curves=(), event=np.array(['A', 'B', 'C', 'D', 'E']),
            rayp=np.array([0.05, 0.06, 0.07, 0.06, 0.06]),
            periods_s=np.array([0.1, 0.5, 1.0, 2.0]),
            vs_km_s=np.array([[1., 2., 3., 4.],
                              [3., 4., 5., 6.],
                              [100., 90., 100., 100.],
                              [100., 100., np.nan, 100.],
                              [np.inf, 100., 100., 100.]]),
            vs0_km_s=np.array([1., 3., 100., 100., 100.]),
            status=np.array([['ok', 'ok', 'ok', 'ok'],
                             ['ok', 'ok', 'ok', 'ok'],
                             ['ok', 'nonpositive_radial', 'ok', 'ok'],
                             ['ok', 'ok', 'ok', 'ok'],
                             ['ok', 'ok', 'ok', 'ok']]),
            vs0_status=('ok',) * 5, reference='seispy-iter',
        )
        for value in vars(self.result).values():
            if isinstance(value, np.ndarray):
                value.setflags(write=False)

    def tearDown(self):
        plt.close('all')

    def test_one_bad_period_excludes_the_event_from_the_entire_curve_and_band(self):
        fig, ax = self.result.plot()
        self.assertIs(fig, ax.figure)
        self.assertEqual(len(ax.lines), 1)  # No horizontal reference line.
        self.assertEqual(ax.lines[0].get_color(), 'red')
        np.testing.assert_array_equal(ax.lines[0].get_xdata(), self.result.periods_s)
        np.testing.assert_allclose(ax.lines[0].get_ydata(), [2., 3., 4., 5.])
        self.assertEqual(len(ax.collections), 1)
        paths = ax.collections[0].get_paths()
        self.assertEqual(len(paths), 1)
        vertices = np.concatenate([path.vertices for path in paths])
        for x, lower, upper in ((0.1, 1., 3.), (0.5, 2., 4.), (1., 3., 5.), (2., 4., 6.)):
            y = vertices[vertices[:, 0] == x, 1]
            self.assertAlmostEqual(y.min(), lower)
            self.assertAlmostEqual(y.max(), upper)

    def test_existing_axes_and_title_are_supported_without_new_figure(self):
        fig, ax = plt.subplots()
        count = len(plt.get_fignums())
        returned_fig, returned_ax = self.result.plot(ax=ax, title='XX.TEST')
        self.assertIs(returned_fig, fig)
        self.assertIs(returned_ax, ax)
        self.assertEqual(len(plt.get_fignums()), count)
        self.assertEqual(ax.get_title(), 'XX.TEST')
        self.assertEqual(ax.get_xlabel(), 'Period (s)')
        self.assertIn('km/s', ax.get_ylabel())
        np.testing.assert_allclose(ax.get_xlim(), [0., 2.])

    def test_no_valid_measurements_raise_before_creating_a_figure(self):
        for values in (np.full((5, 4), np.nan), self.result.vs_km_s):
            with self.subTest(all_nonfinite=np.isnan(values).all()):
                result = replace(self.result, vs_km_s=values,
                                 status=np.full((5, 4), 'invalid_vertical'))
                count = len(plt.get_fignums())
                with self.assertRaisesRegex(VsAppError, 'No event is valid'):
                    result.plot()
                self.assertEqual(len(plt.get_fignums()), count)

    def test_each_period_has_values_but_no_complete_event_is_rejected(self):
        status = self.result.status.copy()
        status[0, 0] = 'invalid_vertical'
        status[1, 3] = 'invalid_vertical'
        result = replace(self.result, status=status)
        # There is still at least one valid value in every column.
        self.assertTrue(np.all(np.any((status == 'ok') & np.isfinite(result.vs_km_s), axis=0)))
        count = len(plt.get_fignums())
        with self.assertRaisesRegex(VsAppError, 'No event is valid'):
            result.plot()
        self.assertEqual(len(plt.get_fignums()), count)

    def test_one_complete_event_has_zero_width_band(self):
        status = self.result.status.copy()
        status[1, 2] = 'invalid_vertical'
        result = replace(self.result, status=status)
        _, ax = result.plot()
        np.testing.assert_array_equal(ax.lines[0].get_ydata(), [1., 2., 3., 4.])
        vertices = ax.collections[0].get_paths()[0].vertices
        for x, expected in zip(result.periods_s, [1., 2., 3., 4.], strict=True):
            np.testing.assert_array_equal(vertices[vertices[:, 0] == x, 1], expected)

    def test_plot_preserves_measurements_and_qc_arrays(self):
        before = {key: value.copy() for key, value in vars(self.result).items()
                  if isinstance(value, np.ndarray)}
        self.result.plot()
        for key, value in before.items():
            np.testing.assert_array_equal(getattr(self.result, key), value)
            self.assertFalse(getattr(self.result, key).flags.writeable)

    def test_display_is_opt_in(self):
        with patch('matplotlib.pyplot.show') as show:
            self.result.plot()
            show.assert_not_called()
            self.result.plot(show=True)
            show.assert_called_once_with()


if __name__ == '__main__':
    unittest.main()
