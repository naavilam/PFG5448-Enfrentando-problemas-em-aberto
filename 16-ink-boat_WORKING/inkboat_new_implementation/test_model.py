import unittest
from dataclasses import replace
import numpy as np
from model import (Params, WATER_VISCOSITY_PA_S, capillary_force,
                   constant_force_response, contact_line, simulate)


class ForceBalanceTests(unittest.TestCase):
    def test_uniform_gamma_cancels_for_different_perimeters_and_resolutions(self):
        for perimeter in (0.002, 0.02, 0.1):
            for gamma in (0.0, 0.030, 0.100):
                for count in (16, 256, 1024):
                    p = Params(perimeter_m=perimeter, gamma_n_m=gamma, contour_points=count)
                    np.testing.assert_allclose(capillary_force(p), [0, 0], atol=1e-18)

    def test_known_first_harmonic_has_correct_vector_force(self):
        # Analytic integral for gamma0 + g1*cos(theta): Fx=P*g1/2, Fy=0.
        p = Params()
        theta, _, _ = contact_line(p)
        amplitude = 0.01
        gamma = p.gamma_n_m + amplitude*np.cos(theta)
        expected = p.perimeter_m*amplitude/2
        np.testing.assert_allclose(capillary_force(p, gamma), [expected, 0], rtol=1e-12, atol=1e-18)
        gamma = p.gamma_n_m + amplitude*np.sin(theta)
        np.testing.assert_allclose(capillary_force(p, gamma), [0, expected], rtol=1e-12, atol=1e-18)

    def test_rest_remains_at_rest_without_external_force(self):
        r = simulate(Params(duration_s=60))
        for key in ('velocity_m_s', 'acceleration_m_s2', 'position_m', 'total_force_n'):
            np.testing.assert_array_equal(r[key], 0)

    def test_loaded_solution_initial_condition_terminal_speed_and_force_balance(self):
        p = Params(external_force_n=1e-8, duration_s=60)
        r = simulate(p)
        self.assertEqual(r['velocity_m_s'][0], 0)
        self.assertEqual(r['position_m'][0], 0)
        self.assertAlmostEqual(r['acceleration_m_s2'][0], p.external_force_n/p.mass_kg)
        self.assertAlmostEqual(r['velocity_m_s'][-1], p.external_force_n/p.drag_kg_s)
        np.testing.assert_allclose(p.mass_kg*r['acceleration_m_s2'], r['total_force_n'], atol=1e-22)
        self.assertTrue((r['velocity_m_s'] >= 0).all())
        self.assertTrue((np.diff(r['velocity_m_s']) >= 0).all())

    def test_signed_force_and_inertial_limit(self):
        time = np.array([0., 0.1, 1., 10.])
        v, a, x = constant_force_response(time, 0.2, 0., -0.4)
        np.testing.assert_allclose(v, -2*time)
        np.testing.assert_allclose(a, -2)
        np.testing.assert_allclose(x, -time**2)
        v, a, x = constant_force_response(time, 0.2, 0.01, -0.4)
        self.assertTrue((v <= 0).all() and (a <= 0).all() and (x <= 0).all())

    def test_small_time_displacement_is_stable(self):
        t = np.array([1e-14, 1e-12, 1e-10])
        _, _, x = constant_force_response(t, 1.0, 1.0, 1.0)
        np.testing.assert_allclose(x, 0.5*t**2, rtol=1e-9, atol=0)

    def test_geometry_and_drag_units(self):
        p = Params(perimeter_m=2*np.pi*0.003, shear_depth_m=0.001)
        self.assertAlmostEqual(p.radius_m, 0.003)
        self.assertAlmostEqual(p.wetted_area_m2, np.pi*0.003**2)
        self.assertAlmostEqual(p.drag_kg_s, WATER_VISCOSITY_PA_S*np.pi*0.003**2/0.001)
        self.assertAlmostEqual(replace(p, perimeter_m=2*p.perimeter_m).drag_kg_s, 4*p.drag_kg_s)

    def test_invalid_parameters_fail(self):
        for kwargs in ({'mass_kg':0}, {'perimeter_m':-1}, {'gamma_n_m':float('nan')},
                       {'shear_depth_m':0}, {'points':1}, {'external_force_n':float('inf')}):
            with self.assertRaises(ValueError):
                Params(**kwargs)


if __name__ == '__main__':
    unittest.main()
