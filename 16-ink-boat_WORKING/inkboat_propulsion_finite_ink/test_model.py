import unittest
from dataclasses import replace
import numpy as np
from model import (Params, WATER_VISCOSITY_PA_S, capillary_force,
                   constant_force_response, contact_line, simulate, ink_fraction, finite_ink_response, ink_active)


class ForceBalanceTests(unittest.TestCase):
    def test_uniform_gamma_cancels_for_different_perimeters_and_resolutions(self):
        for perimeter in (0.002, 0.02, 0.1):
            for gamma in (0.0, 0.030, 0.100):
                for count in (16, 256, 1024):
                    p = Params(perimeter_m=perimeter, gamma_n_m=gamma, gamma_contrast_n_m=0, contour_points=count)
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

    def test_rest_remains_at_rest_without_contrast(self):
        r = simulate(Params(duration_s=60, gamma_contrast_n_m=0))
        for key in ('velocity_m_s', 'acceleration_m_s2', 'position_m', 'total_force_n'):
            np.testing.assert_array_equal(r[key], 0)

    def test_loaded_solution_initial_condition_terminal_speed_and_force_balance(self):
        p = Params(duration_s=60, gamma_contrast_n_m=2e-6, ejection_rate_m3_s=0)
        force = p.perimeter_m*p.gamma_contrast_n_m/4
        r = simulate(p)
        self.assertEqual(r['velocity_m_s'][0], 0)
        self.assertEqual(r['position_m'][0], 0)
        self.assertAlmostEqual(r['acceleration_m_s2'][0], force/p.mass_kg)
        self.assertAlmostEqual(r['velocity_m_s'][-1], force/p.drag_kg_s)
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

    def test_prescribed_contrast_produces_correct_force_and_motion(self):
        for contrast in (1e-5, -1e-5, 0.0):
            p = Params(gamma_contrast_n_m=contrast)
            expected = p.perimeter_m*contrast/4
            np.testing.assert_allclose(capillary_force(p), [expected, 0], atol=1e-18)
            r = simulate(p)
            self.assertAlmostEqual(r['acceleration_m_s2'][0], expected/p.mass_kg)
            np.testing.assert_allclose(p.mass_kg*r['acceleration_m_s2'], r['total_force_n'], atol=1e-21)
            if contrast:
                self.assertGreater(r['velocity_m_s'][-1]*contrast, 0)
        a = capillary_force(Params(gamma_n_m=0.03))
        b = capillary_force(Params(gamma_n_m=0.09))
        np.testing.assert_allclose(a, b, atol=1e-18)
        with self.assertRaises(ValueError):
            Params(gamma_n_m=0.01, gamma_contrast_n_m=0.03)

    def test_ink_budget_and_contrast_decay(self):
        p = Params()
        time = np.array([0, p.empty_time_s/2, p.empty_time_s, 2*p.empty_time_s])
        np.testing.assert_allclose(ink_fraction(time, p), [1, 0.5, 0, 0])
        for t, fraction in zip(time, [1, 1, 0, 0]):
            np.testing.assert_allclose(capillary_force(p, time_s=t),
                                       [p.perimeter_m*p.gamma_contrast_n_m*fraction/4, 0], atol=1e-18)
        r = simulate(replace(p, duration_s=10))
        np.testing.assert_allclose(r['ink_remaining_m3']+r['ink_ejected_m3'], p.ink_volume_m3, atol=1e-26)
        self.assertTrue((r['ink_remaining_m3'] >= 0).all())
        self.assertTrue((np.diff(r['ink_remaining_m3']) <= 0).all())
        self.assertIn(p.empty_time_s, r['time_s'])
        self.assertEqual(replace(p, ejection_rate_m3_s=2*p.ejection_rate_m3_s).empty_time_s, p.empty_time_s/2)

    def test_no_volume_and_no_consumption_limits(self):
        p = Params(ink_volume_m3=0)
        r = simulate(p)
        for key in ('contrast_n_m', 'capillary_force_n', 'velocity_m_s', 'acceleration_m_s2'):
            np.testing.assert_array_equal(r[key], 0)
        p = Params(ejection_rate_m3_s=0)
        r = simulate(p)
        expected = constant_force_response(r['time_s'], p.mass_kg, p.drag_kg_s,
                                           p.perimeter_m*p.gamma_contrast_n_m/4)
        for key, values in zip(('velocity_m_s','acceleration_m_s2','position_m'), expected):
            np.testing.assert_allclose(r[key], values)
        np.testing.assert_array_equal(r['ink_remaining_m3'], p.ink_volume_m3)

    def test_continuity_at_depletion_and_viscous_coasting(self):
        p = Params()
        T, tau = p.empty_time_s, p.relaxation_time_s
        time = np.array([T-1e-8, T, T+1e-8, T+tau])
        v, a, x = finite_ink_response(time, p)
        np.testing.assert_allclose(v[:3], v[1], rtol=1e-7)
        self.assertAlmostEqual(a[1]-a[0], -p.perimeter_m*p.gamma_contrast_n_m/(4*p.mass_kg), places=9)
        np.testing.assert_allclose(x[:3], x[1], rtol=1e-7)
        self.assertGreater(v[1], 0)
        self.assertLess(a[1], 0)
        self.assertAlmostEqual(v[-1]/v[1], np.exp(-1))
        # Negative contrast reverses the entire trajectory.
        nv, na, nx = finite_ink_response(time, replace(p, gamma_contrast_n_m=-p.gamma_contrast_n_m))
        np.testing.assert_allclose(nv, -v)
        np.testing.assert_allclose(na, -a)
        np.testing.assert_allclose(nx, -x)

    def test_cutoff_solution_matches_independent_rk4_integration(self):
        p = Params()
        total = 2*p.empty_time_s
        steps = 4000
        dt = total/steps
        state = np.zeros(2)  # velocity and displacement
        force0 = p.perimeter_m*p.gamma_contrast_n_m/4
        # Integrate separate constant-force intervals so RK4 never straddles
        # the discontinuity or averages the two one-sided accelerations.
        def rhs(y, force):
            return np.array([(force-p.drag_kg_s*y[0])/p.mass_kg, y[0]])
        for force in (force0, 0):
            for k in range(steps//2):
                k1 = rhs(state, force)
                k2 = rhs(state+dt*k1/2, force)
                k3 = rhs(state+dt*k2/2, force)
                k4 = rhs(state+dt*k3, force)
                state += dt*(k1+2*k2+2*k3+k4)/6
        v, a, x = finite_ink_response(np.array([total]), p)
        np.testing.assert_allclose([v[0], x[0]], state, rtol=1e-10, atol=1e-16)

    def test_fast_depletion_and_tiny_times_are_stable(self):
        p = Params(ejection_rate_m3_s=1e-3)
        time = np.array([0, p.empty_time_s/2, p.empty_time_s, 1e-6])
        v, a, x = finite_ink_response(time, p)
        self.assertTrue(np.isfinite(v).all() and np.isfinite(x).all())
        self.assertTrue((v >= 0).all() and (x >= 0).all())

    def test_invalid_parameters_fail(self):
        for kwargs in ({'mass_kg':0}, {'perimeter_m':-1}, {'gamma_n_m':float('nan')},
                       {'shear_depth_m':0}, {'points':1}, {'gamma_contrast_n_m':float('inf')}, {'ink_volume_m3':-1}, {'ejection_rate_m3_s':-1}):
            with self.assertRaises(ValueError):
                Params(**kwargs)


if __name__ == '__main__':
    unittest.main()
