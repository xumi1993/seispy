"""INI parameters, backward-compatible keyword calls and SeisPy logging."""

import logging
import tempfile
import unittest
from dataclasses import fields
from pathlib import Path
from unittest.mock import patch

import numpy as np
import test_vsapp_inversion as fixtures
from seispy.inversion import InvPara, invert_station_vsapp, invert_vsapp, invpara
from seispy.inversion.forward import VsappForward
from seispy.inversion.invpara import resolve_parameters
from seispy.setuplog import SetupLog, inversion_logger
from test_vsapp_inversion import FORWARD, PERIODS, make_model, native_curve


class TestInvPara(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / 'inversion.cfg'

    def write_config(self, content):
        self.path.write_text(content)
        return self.path

    def test_sectioned_file_reads_types_arrays_and_all_parameters(self):
        example = Path(__file__).parents[1] / 'seispy/data/vsapp_inversion.cfg'
        parameters = invpara(example)
        np.testing.assert_allclose(parameters.periods_s, PERIODS)
        self.assertEqual(parameters.rayp, .06)
        self.assertEqual(parameters.f0, 3.5)
        self.assertIs(parameters.zero_halfspace, True)
        self.assertEqual(parameters.misfit_window, 5)
        self.assertEqual(set(parameters.kernel_kwargs),
                         {'dt', 'npts', 'shift', 'pre_filt', 'itmax', 'minderr'})
        # The distributed example documents every setting in the parameter class.
        text = example.read_text()
        for item in fields(InvPara):
            self.assertIn(item.name + ' =', text)
        self.assertIn('optimizer: gd', str(parameters))
        path = self.write_config('[data]\nsigma=.03 .04\n'
                                 '[forward]\nmethod=water\npre_filt=.1, 2\n'
                                 '[inversion]\nzero_halfspace=no\n')
        parameters = InvPara.read_para(path)
        np.testing.assert_array_equal(parameters.sigma, [.03, .04])
        self.assertFalse(parameters.zero_halfspace)
        self.assertIn('wlevel', parameters.kernel_kwargs)
        self.assertNotIn('itmax', parameters.kernel_kwargs)

    def test_scalar_uncertainty_and_blank_defaults(self):
        parameters = invpara(self.write_config('[data]\nsigma=.1\n[forward]\nf0=\n'))
        self.assertEqual(parameters.sigma, .1)
        self.assertIsNone(parameters.f0)
        self.assertEqual(parameters.max_step_lnvs, InvPara().max_step_lnvs)

    def test_unknown_or_malformed_configuration_is_rejected(self):
        cases = ['[wrong]\nf0=3.5', '[inversion]\nmisfit_toll=.01',
                 '[forward]\nnpts=3.2', '[forward]\nmethod=bad',
                 '[data]\nsigma=__import__("os")', '[inversion',
                 '[inversion]\nmisfit_window=0']
        for content in cases:
            with self.subTest(content=content):
                with self.assertRaises((ValueError, TypeError)):
                    invpara(self.write_config(content))
        with self.assertRaises(FileNotFoundError):
            invpara(Path(self.directory.name) / 'missing.cfg')

    def test_parameter_copy_and_override_precedence(self):
        original = InvPara(f0=3.5, rayp=.06, periods_s=PERIODS, sigma=np.ones(7))
        result = resolve_parameters(original, dt=.1, kernel_kwargs={'dt': .2, 'npts': 512})
        self.assertEqual(result.dt, .1)
        self.assertEqual(result.npts, 512)
        self.assertEqual(original.dt, .05)
        self.assertEqual(original.npts, 1024)
        result.sigma[0] = 9.
        self.assertEqual(original.sigma[0], 1.)
        with self.assertRaises(ValueError):
            resolve_parameters(original, self.path)
        with self.assertRaises(TypeError):
            resolve_parameters(original, unknown_option=1)

    def test_file_object_and_keyword_calls_have_identical_results(self):
        self.write_config('[data]\nperiods_s=' + ','.join(map(str, PERIODS)) + '\n'
                          '[forward]\nf0=3.5\nrayp=.06\nmethod=water\n'
                          'dt=.05\nnpts=512\nshift=5\n'
                          '[inversion]\nmax_iterations=1\n')
        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'water')
        parameters = invpara(self.path)
        by_file = invert_vsapp(make_model(), observed_vsapp=observed, cfg_file=self.path)
        by_object = invert_vsapp(make_model(), observed_vsapp=observed, para=parameters)
        by_keywords = invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5,
                                    method='water', max_iterations=1, kernel_kwargs=FORWARD)
        for other in (by_object, by_keywords):
            np.testing.assert_array_equal(other.model.vs, by_file.model.vs)
            np.testing.assert_array_equal(other.misfit_history, by_file.misfit_history)
            self.assertEqual(other.n_evaluations, by_file.n_evaluations)
        unchanged = parameters.periods_s.copy()
        override = invert_vsapp(make_model(), PERIODS, observed, cfg_file=self.path,
                                 max_iterations=0)
        self.assertEqual(override.n_iterations, 0)
        np.testing.assert_array_equal(parameters.periods_s, unchanged)
        self.assertEqual(parameters.max_iterations, 1)

    def test_station_uses_measurement_periods_and_mean_rayp_without_mutating_parameters(self):
        fixture = fixtures.TestStationInversion()
        fixture.setUp()
        parameters = InvPara(f0=3.5, rayp=.09, periods_s=np.array([1., 2.]),
                             max_iterations=0, **FORWARD)
        result = invert_station_vsapp(make_model(), fixture.measurements, para=parameters)
        self.assertAlmostEqual(result.rayp, .06)
        np.testing.assert_array_equal(result.periods_s, PERIODS)
        self.assertEqual(parameters.rayp, .09)
        np.testing.assert_array_equal(parameters.periods_s, [1., 2.])

    def test_invalid_inputs_fail_before_forward_and_are_logged(self):
        with self.assertLogs('Inv', level='ERROR') as messages:
            with patch.object(VsappForward, 'evaluate') as forward:
                with self.assertRaises(ValueError):
                    invert_vsapp(make_model(), PERIODS, np.ones(7), f0=3.5)
                forward.assert_not_called()
        self.assertIn('rayp', '\n'.join(messages.output))


class TestInversionLogging(unittest.TestCase):
    def test_stream_setup_does_not_open_files_or_duplicate_handlers(self):
        logger = logging.getLogger('Inv')
        with patch('seispy.setuplog.logging.FileHandler') as handler:
            setup = SetupLog(filename=None)
            before = tuple(logger.handlers)
            self.assertIs(inversion_logger(setup), logger)
            self.assertIs(inversion_logger(), logger)
            self.assertIs(inversion_logger(), logger)
            handler.assert_not_called()
        self.assertEqual(tuple(logger.handlers), before)

    def test_injected_logger_receives_progress_trials_and_stop_reason_without_print(self):
        logger = logging.getLogger('test.vsapp.custom')
        observed = native_curve(make_model((1.2, 2.5, 3.5)), 'water')
        with self.assertLogs(logger, level='DEBUG') as messages:
            with patch('builtins.print') as print_mock:
                result = invert_vsapp(make_model(), PERIODS, observed, rayp=.06, f0=3.5,
                                      method='water', max_iterations=1,
                                      kernel_kwargs=FORWARD, log=logger)
                print_mock.assert_not_called()
        output = '\n'.join(messages.output)
        for text in ('Starting Vsapp inversion', 'Forward/kernel evaluation', 'Iteration 0:',
                     'Iteration 1:', 'misfit=', 'step_length=', result.message):
            self.assertIn(text, output)
        self.assertTrue(any(record.levelno == logging.WARNING for record in messages.records))

    def test_convergence_and_station_qc_are_logged(self):
        fixture = fixtures.TestStationInversion()
        fixture.setUp()
        with self.assertLogs('Inv', level='INFO') as messages:
            result = invert_station_vsapp(make_model(), fixture.measurements,
                                           f0=3.5, kernel_kwargs=FORWARD)
        self.assertTrue(result.converged)
        output = '\n'.join(messages.output)
        self.assertIn('retained 2/3 events', output)
        self.assertIn('Log-gradient tolerance reached', output)

    def test_explicit_log_level_is_preserved(self):
        logger = logging.getLogger('Inv')
        before = logger.level
        try:
            logger.setLevel(logging.DEBUG)
            self.assertIs(inversion_logger(), logger)
            self.assertEqual(logger.level, logging.DEBUG)
        finally:
            logger.setLevel(before)
