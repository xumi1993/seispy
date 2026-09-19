"""RFStation adapters use SeisPy metadata without changing the stored RFs."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from obspy import Stream, Trace
from obspy.io.sac import SACTrace
from seispy import DepModel, RFStation, SynSeis
from seispy.decon import RFTrace
from seispy.geo import skm2srad
from seispy.vsapp import (
    VsAppError,
    gaussian_vertical,
    seispy_iter_vertical,
)
from seispy.vsapp import (
    _compute_vsapp_from_components as compute_paired,
)


class TestStationVsApp(unittest.TestCase):
    def native_iterative_pair(self, dt, shift, npts=512):
        waveform = np.exp(-((np.arange(npts) * dt - 3.0) / 0.4) ** 2)
        z = Trace(waveform, header={'delta': dt})
        r = Trace(0.3 * waveform, header={'delta': dt})
        return (
            RFTrace.deconvolve(r, z, method='iter', tshift=shift, f0=2.0),
            RFTrace.deconvolve(z, z, method='iter', tshift=shift, f0=2.0),
        )

    def make_station(self, normalization='area'):
        time = np.arange(301) * 0.05 - 5.0
        rayp = np.array([0.04, 0.08])
        ratios = np.array([0.3, 0.5])
        factors = np.array([2.0, 3.5])
        vertical = np.stack([
            gaussian_vertical(time, f0, normalization=normalization)
            for f0 in factors
        ])
        stream = Stream([
            Trace(ratio * z, header={'delta': 0.05, 'tshift': 5.0, 'f0': f0})
            for ratio, z, f0 in zip(ratios, vertical, factors, strict=True)
        ])
        station = RFStation.read_stream(stream, rayp, np.array([200.0, 10.0]))
        station.event[:] = ['A', 'B']
        station.phase[:] = 'P'
        return station, vertical, rayp, ratios

    def test_each_event_uses_its_own_slowness_and_gaussian_factor(self):
        station, vertical, rayp, ratios = self.make_station()
        result = station.compute_vsapp([0.1, 0.7, 3.0], reference='gaussian-area')
        expected = np.sin(np.arctan(ratios) / 2) / rayp
        np.testing.assert_allclose(result.vs_km_s, np.repeat(expected[:, None], 3, axis=1))
        np.testing.assert_allclose(result.vs0_km_s, expected)
        np.testing.assert_allclose(result.rayp, rayp)
        np.testing.assert_allclose(station.rayp, skm2srad(rayp))
        np.testing.assert_array_equal(result.event, ['A', 'B'])
        self.assertTrue(np.all(result.status == 'ok'))
        self.assertEqual(result.vs0_status, ('ok', 'ok'))
        for i, curve in enumerate(result.curves):
            direct = compute_paired(station.time_axis, station.data_prime[i], vertical[i],
                                   rayp=rayp[i], periods_s=[0.1, 0.7, 3.0])
            np.testing.assert_allclose(curve.radial_integral, direct.radial_integral)

    def test_public_station_api_has_no_vertical_input_and_requires_f0(self):
        station, vertical, _, _ = self.make_station()
        with self.assertRaises(TypeError):
            station.compute_vsapp([0.2, 1.0], vertical=vertical)
        station.f0[:] = np.nan
        with self.assertRaises(VsAppError):
            station.compute_vsapp([0.2, 1.0])

    def test_invalid_velocity_keeps_event_and_qc_status(self):
        station, vertical, _, _ = self.make_station()
        station.data_prime[1] *= -1
        result = station.compute_vsapp([0.2, 1.0], reference='gaussian-area')
        self.assertTrue(np.isfinite(result.vs_km_s[0]).all())
        self.assertTrue(np.isnan(result.vs_km_s[1]).all())
        self.assertTrue(np.all(result.status[1] == 'nonpositive_radial'))
        self.assertEqual(result.vs0_status[1], 'nonpositive_radial')
        np.testing.assert_array_equal(result.event, ['A', 'B'])

    def test_peak_reference_uses_matched_amplitude_convention(self):
        station, _, rayp, ratios = self.make_station(normalization='peak')
        result = station.compute_vsapp([0.2, 1.0], reference='gaussian-peak')
        expected = np.sin(np.arctan(ratios) / 2) / rayp
        np.testing.assert_allclose(result.vs0_km_s, expected)
        np.testing.assert_allclose(result.vs_km_s[:, 0], expected)

    def test_sort_order_and_snapshot_ownership(self):
        station, _, _, _ = self.make_station()
        before_data = station.data_prime.copy()
        result = station.compute_vsapp([0.2, 1.0], reference='gaussian-area')
        np.testing.assert_array_equal(station.data_prime, before_data)
        station.sort('bazi')
        sorted_result = station.compute_vsapp([0.2, 1.0], reference='gaussian-area')
        np.testing.assert_array_equal(sorted_result.event, ['B', 'A'])
        np.testing.assert_allclose(sorted_result.vs_km_s, result.vs_km_s[::-1])
        np.testing.assert_array_equal(result.event, ['A', 'B'])
        for array in (result.event, result.rayp, result.periods_s, result.vs_km_s,
                      result.vs0_km_s, result.status):
            self.assertFalse(array.flags.writeable)

    def test_native_synthetic_stream_matches_actual_vertical_deconvolution(self):
        model = DepModel.read_layer_model(
            np.array([0.0, 1.0, 5.0]), [1.0, 4.0, 0.0],
            [2.0, 4.5, 6.3], [0.8, 2.5, 3.7], [2.0, 2.5, 2.8],
        )
        rayp = np.array([0.045, 0.075])
        synth = SynSeis(model, rayp, 0.05, npts=512)
        synth.run_fwd()
        stream = synth.run_deconvolution(pre_filt=None, shift=5.0, f0=2.0, method='iter')
        vertical = np.stack([
            RFTrace.deconvolve(z, z, method='iter', tshift=5.0, f0=2.0).data
            for z in synth.zstream
        ])
        station = RFStation.read_stream(stream, rayp, np.array([10.0, 20.0]))
        reconstructed = station.compute_vsapp([0.1, 0.5, 2.0])
        self.assertEqual(reconstructed.reference, 'seispy-iter')
        for i, p in enumerate(rayp):
            direct = compute_paired(station.time_axis, stream[i].data, vertical[i],
                                   rayp=p, periods_s=[0.1, 0.5, 2.0])
            np.testing.assert_allclose(reconstructed.vs_km_s[i], direct.vs_km_s,
                                       rtol=2e-12, atol=2e-12)

    def test_truncated_iterative_rf_accepts_original_deconvolution_length(self):
        full_time = np.arange(1024) * 0.05 - 5.0
        original_z = seispy_iter_vertical(full_time, 2.0)
        saved_z = original_z[:301]
        stream = Stream([Trace(0.3 * saved_z, header={
            'delta': 0.05, 'tshift': 5.0, 'f0': 2.0,
        })])
        station = RFStation.read_stream(stream, 0.06, 20.0)
        result = station.compute_vsapp([0.1, 0.7, 2.0], reference='seispy-iter',
                                       deconvolution_npts=1024)
        direct = compute_paired(station.time_axis, station.data_prime[0], saved_z,
                                rayp=0.06, periods_s=[0.1, 0.7, 2.0])
        np.testing.assert_allclose(result.vs_km_s[0], direct.vs_km_s)
        with self.assertRaises(VsAppError):
            station.compute_vsapp([0.1], reference='seispy-iter', deconvolution_npts=300)

    def test_custom_rftrace_without_numeric_iteration_remains_readable(self):
        radial, _ = self.native_iterative_pair(0.05, 5.0)
        for value in (None, 'unavailable'):
            with self.subTest(iteration=value):
                trace = radial.copy()
                trace.stats.iter = value
                station = RFStation.read_stream(Stream([trace]), 0.06, 20.0)
                self.assertFalse(hasattr(station, '_vsapp_iter_timing'))
                np.testing.assert_array_equal(station.data_prime[0], radial.data)

    def test_native_shift_truncation_boundary_retains_original_index(self):
        radial, vertical = self.native_iterative_pair(0.1, 0.3)
        station = RFStation.read_stream(Stream([radial]), 0.06, 20.0)
        self.assertEqual(station._vsapp_iter_timing, (0.1, 0.3, 2))
        before = station.data_prime.copy()
        actual = compute_paired(station.time_axis, radial.data, vertical.data,
                                rayp=0.06, periods_s=[0.22, 0.25])
        reference = station.compute_vsapp([0.22, 0.25])
        np.testing.assert_allclose(reference.vs_km_s[0], actual.vs_km_s, rtol=2e-12, atol=2e-12)
        nominal = station.compute_vsapp([0.22, 0.25], reference='seispy-iter',
                                        deconvolution_shift_samples=3)
        self.assertGreater(np.max(np.abs(nominal.vs_km_s - actual.vs_km_s)), 1e-3)
        np.testing.assert_array_equal(station.data_prime, before)

    def test_explicit_shift_index_accepts_non_grid_arrival(self):
        radial, vertical = self.native_iterative_pair(0.1, 0.35)
        times = np.arange(len(radial)) * 0.1 - 0.35
        rounded = int(round(-times[0] / np.median(np.diff(times))))
        np.testing.assert_array_equal(
            seispy_iter_vertical(times, 2.0),
            seispy_iter_vertical(times, 2.0, deconvolution_shift_samples=rounded),
        )
        reconstructed = seispy_iter_vertical(times, 2.0, deconvolution_shift_samples=3)
        np.testing.assert_allclose(reconstructed, vertical.data, rtol=2e-12, atol=2e-12)
        for invalid in (-1, 512, True, np.bool_(False), 1.5):
            with self.subTest(index=invalid), self.assertRaises(VsAppError):
                seispy_iter_vertical(times, 2.0, deconvolution_shift_samples=invalid)

    def test_sac_station_automatically_rounds_non_grid_arrivals(self):
        for shift, index in ((5.00002, 100), (5.019, 100), (5.031, 101)):
            with self.subTest(shift=shift), TemporaryDirectory() as directory:
                path = Path(directory)
                radial, _ = self.native_iterative_pair(0.05, shift)
                SACTrace(
                    data=radial.data.astype(np.float32), delta=0.05, b=-shift,
                    knetwk='XX', kstnm='TEST', stla=0.0, stlo=0.0, stel=0.0,
                    user0=0.06, user1=2.0,
                ).write(str(path / 'EVENT_P_R.sac'))
                (path / 'XX.TESTfinallist.dat').write_text(
                    'EVENT P 0 0 10 60 20 0.06 6 2.0\n'
                )
                station = RFStation(str(path), only_r=True)
                self.assertFalse(hasattr(station, '_vsapp_iter_timing'))
                before_time = station.time_axis.copy()
                before_data = station.data_prime.copy()
                actual = station.compute_vsapp([0.11, 0.7, 2.0])
                expected = station.compute_vsapp(
                    [0.11, 0.7, 2.0], deconvolution_shift_samples=index,
                )
                np.testing.assert_array_equal(actual.vs_km_s, expected.vs_km_s)
                np.testing.assert_array_equal(actual.vs0_km_s, expected.vs0_km_s)
                self.assertTrue(np.all(actual.status == 'ok'))
                np.testing.assert_array_equal(station.time_axis, before_time)
                np.testing.assert_array_equal(station.data_prime, before_data)

    def test_sac_station_nominal_and_explicit_original_shift(self):
        for dt, shift, periods, index in (
            (0.05, 5.0, [0.11, 1.0, 3.0], None),
            (0.1, 0.3, [0.22, 0.25], 2),
        ):
            with self.subTest(dt=dt, shift=shift), TemporaryDirectory() as directory:
                path = Path(directory)
                radial, vertical = self.native_iterative_pair(dt, shift)
                for component, trace in (('R', radial), ('Z', vertical)):
                    SACTrace(
                        data=trace.data.astype(np.float32), delta=dt, b=-shift,
                        knetwk='XX', kstnm='TEST', stla=0.0, stlo=0.0, stel=0.0,
                        user0=0.06, user1=2.0,
                    ).write(str(path / f'EVENT_P_{component}.sac'))
                (path / 'XX.TESTfinallist.dat').write_text(
                    'EVENT P 0 0 10 60 20 0.06 6 2.0\n'
                )
                station = RFStation(str(path), only_r=True)
                self.assertFalse(hasattr(station, '_vsapp_iter_timing'))
                measured_z = SACTrace.read(str(path / 'EVENT_P_Z.sac')).data
                measured = compute_paired(station.time_axis, station.data_prime[0], measured_z,
                                          rayp=0.06, periods_s=periods)
                assumed = station.compute_vsapp(periods, reference='seispy-iter',
                                                deconvolution_shift_samples=index)
                np.testing.assert_allclose(assumed.vs_km_s[0], measured.vs_km_s,
                                           rtol=3e-6, atol=3e-6)
                np.testing.assert_allclose(assumed.rayp, [0.06], rtol=1e-6)
                self.assertEqual(assumed.event[0], 'EVENT')

    def test_unknown_references_and_conflicting_options_are_rejected(self):
        station, _, _, _ = self.make_station()
        cases = [
            {'reference': None}, {'reference': 'measured'}, {'reference': 'unknown'},
            {'reference': 'gaussian-area', 'deconvolution_npts': 1024},
            {'reference': 'gaussian-peak', 'deconvolution_shift_samples': 100},
            {'deconvolution_shift_samples': -1},
        ]
        before = station.data_prime.copy()
        for options in cases:
            with self.subTest(options=options), self.assertRaises(VsAppError):
                station.compute_vsapp([0.2, 1.0], **options)
        np.testing.assert_array_equal(station.data_prime, before)
        self.assertFalse(hasattr(station, 'dataz'))

    def test_non_radial_and_s_phase_are_rejected(self):
        station, vertical, _, _ = self.make_station()
        for comp in ('Q', 'L', 'Z', 'T'):
            station.comp = comp
            with self.subTest(comp=comp), self.assertRaises(VsAppError):
                station.compute_vsapp([0.2])
        station.comp = 'R'
        station.phase[1] = 'S'
        with self.assertRaises(VsAppError):
            station.compute_vsapp([0.2])
        station.phase[:] = 'P'
        station.prime_phase = 'S'
        with self.assertRaises(VsAppError):
            station.compute_vsapp([0.2])

    def test_invalid_metadata_and_windows_are_rejected_without_waveform_mutation(self):
        for attr, value in (
            ('rayp', np.array([0.0, 1.0])),
            ('rayp', np.array([1.0])),
            ('rayp', np.array(['bad', 'data'])),
            ('f0', np.array([0.0, 2.0])),
            ('f0', np.array([2.0])),
            ('phase', np.array(['P'])),
            ('event', np.array(['A'])),
            ('ev_num', 3),
        ):
            station, _, _, _ = self.make_station()
            before = station.data_prime.copy()
            setattr(station, attr, value)
            with self.subTest(attr=attr, value=value), self.assertRaises(VsAppError):
                station.compute_vsapp([0.2], reference='gaussian-area')
            np.testing.assert_array_equal(station.data_prime, before)
        station, vertical, _, _ = self.make_station()
        for periods in ([0.01], [6.0], [1.0, 0.5]):
            with self.subTest(periods=periods), self.assertRaises(VsAppError):
                station.compute_vsapp(periods)
        with self.assertRaises(VsAppError):
            station.compute_vsapp([0.2], denominator_rtol=1.0)


if __name__ == '__main__':
    unittest.main()
