"""Equation 3.5 with a finite-ink, time-dependent surface-tension contrast and an explicit thin-disk shear closure.

SI units internally. This is a reduced force balance, not a fluid-field solver.
See README.md for assumptions, the derivation, and sources.
"""
from dataclasses import dataclass
import numpy as np

# Pure liquid water at 20 deg C, approximately 0.1 MPa.
# IAPWS SR6-08(2011), section 3 / Eq. (7), reference value rounded to 5 digits.
WATER_TEMPERATURE_C = 20.0
WATER_PRESSURE_PA = 100_000.0
WATER_VISCOSITY_PA_S = 1.0016e-3
WATER_DENSITY_KG_M3 = 998.2
WATER_SOURCE = "https://iapws.org/technical-guidance/release/LiquidWater"


@dataclass(frozen=True)
class Params:
    mass_kg: float = 0.2e-3
    perimeter_m: float = 20e-3
    gamma_n_m: float = 72e-3  # prescribed constant, NOT an inferred water property
    gamma_contrast_n_m: float = 0.01e-3  # initial front-minus-rear contrast
    ink_volume_m3: float = 0.02e-9  # 0.02 microlitres
    ejection_rate_m3_s: float = 0.005e-9  # 0.005 microlitres/second
    shear_depth_m: float = 0.5e-3
    duration_s: float = 5.0
    points: int = 501
    contour_points: int = 256

    def __post_init__(self):
        for name in ("mass_kg", "perimeter_m", "shear_depth_m", "duration_s"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not np.isfinite(self.gamma_n_m) or self.gamma_n_m < 0:
            raise ValueError("gamma_n_m must be finite and nonnegative")
        if (not np.isfinite(self.gamma_contrast_n_m)
                or abs(self.gamma_contrast_n_m) > 2*self.gamma_n_m):
            raise ValueError("Contrast must satisfy abs(delta_gamma) <= 2*gamma0")
        for name in ("ink_volume_m3", "ejection_rate_m3_s"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        for name, minimum in (("points", 2), ("contour_points", 16)):
            value = getattr(self, name)
            if not isinstance(value, (int, np.integer)) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")

    @property
    def empty_time_s(self):
        if self.ink_volume_m3 == 0:
            return 0.0
        if self.ejection_rate_m3_s == 0:
            return np.inf
        return self.ink_volume_m3 / self.ejection_rate_m3_s

    @property
    def radius_m(self):
        return self.perimeter_m / (2 * np.pi)

    @property
    def wetted_area_m2(self):
        # Wetted horizontal lower face; thin disk, edge area neglected.
        return np.pi * self.radius_m**2

    @property
    def drag_kg_s(self):
        return WATER_VISCOSITY_PA_S * self.wetted_area_m2 / self.shear_depth_m

    @property
    def relaxation_time_s(self):
        return self.mass_kg / self.drag_kg_s


def contact_line(p):
    """Midpoint samples of the circular contact line and its outward conormal.

    The thesis t is tangent to the water-air interface but NORMAL to the line.
    It is the planar radial unit vector for this geometry.
    """
    theta = (np.arange(p.contour_points) + 0.5) * 2 * np.pi / p.contour_points
    conormal = np.column_stack((np.cos(theta), np.sin(theta)))
    return theta, p.radius_m * conormal, conormal


def gamma_on_contour(theta, time_s, p):
    """Single constitutive hook for a later gamma(theta, t) model.

    Contrast is constant while ink remains, then switches off at depletion.
    Front is +x; lower rear tension mimics the effect of ink.
    """
    return p.gamma_n_m + 0.5*p.gamma_contrast_n_m*ink_active(time_s, p)*np.cos(theta)


def capillary_force(p, gamma_values=None, time_s=0.0):
    """Closed-circle integral F = integral gamma*t dl (N).

    Removing a constant baseline is exact for a closed planar circle because
    integral t dl = 0. It prevents roundoff from manufacturing propulsion.
    The prescribed cosine profile gives Fx=P*delta_gamma/4, Fy=0.
    Optional arrays allow independent mathematical validation of the integral.
    """
    theta, _, conormal = contact_line(p)
    gamma = gamma_on_contour(theta, time_s, p) if gamma_values is None else np.asarray(gamma_values, dtype=float)
    if gamma.shape != theta.shape or not np.isfinite(gamma).all() or (gamma < 0).any():
        raise ValueError("Expected one finite nonnegative gamma per contour point")
    dl = p.perimeter_m / p.contour_points
    return dl * np.sum((gamma - gamma[0])[:, None] * conormal, axis=0)


def constant_force_response(time_s, mass_kg, drag_kg_s, force_n):
    """Exact U(0)=0 solution of m*dU/dt = F - b*U, including b=0 limit.

    Return signed velocity, acceleration and displacement in SI units.
    """
    time = np.asarray(time_s, dtype=float)
    if (not np.isfinite(time).all() or (time < 0).any()
            or not np.isfinite(mass_kg) or mass_kg <= 0
            or not np.isfinite(drag_kg_s) or drag_kg_s < 0
            or not np.isfinite(force_n)):
        raise ValueError("Invalid time, mass, drag or force")
    if drag_kg_s == 0:
        return force_n * time / mass_kg, np.full_like(time, force_n / mass_kg), 0.5 * force_n * time**2 / mass_kg
    tau = mass_kg / drag_kg_s
    q = time / tau
    rise = -np.expm1(-q)
    terminal = force_n / drag_kg_s
    velocity = terminal * rise
    acceleration = (force_n / mass_kg) * np.exp(-q)
    # Stable q + expm1(-q) for very short times.
    factor = q + np.expm1(-q)
    small = q < 1e-4
    factor = np.where(small, q**2 * (0.5 - q/6 + q**2/24 - q**3/120), factor)
    displacement = terminal * tau * factor
    return velocity, acceleration, displacement


def ink_fraction(time_s, p):
    """Remaining volume fraction for constant ejection rate.

    Q=0 keeps the initial contrast if V0>0; V0=0 gives no contrast at all.
    """
    time = np.asarray(time_s, dtype=float)
    if not np.isfinite(time).all() or (time < 0).any():
        raise ValueError("Time must be finite and nonnegative")
    if p.ink_volume_m3 == 0:
        return np.zeros_like(time)
    return np.clip(1 - p.ejection_rate_m3_s*time/p.ink_volume_m3, 0, 1)


def ink_active(time_s, p):
    """Right-continuous cutoff: 1 before V0/Q, 0 at and after depletion."""
    time = np.asarray(time_s, dtype=float)
    if not np.isfinite(time).all() or (time < 0).any():
        raise ValueError("Time must be finite and nonnegative")
    return ((p.ink_volume_m3 > 0) & (time < p.empty_time_s)).astype(float)


def finite_ink_response(time_s, p):
    """Exact response to a constant capillary force cut off at T=V0/Q.

    U and x remain continuous; acceleration jumps by -F0/m at depletion.
    Q=0 gives constant forcing indefinitely when V0>0.
    """
    time = np.asarray(time_s, dtype=float)
    if not np.isfinite(time).all() or (time < 0).any():
        raise ValueError("Time must be finite and nonnegative")
    force0 = p.perimeter_m*p.gamma_contrast_n_m/4
    if p.ink_volume_m3 == 0:
        return np.zeros_like(time), np.zeros_like(time), np.zeros_like(time)
    active_time = np.minimum(time, p.empty_time_s)
    velocity_active, _, position_active = constant_force_response(
        active_time, p.mass_kg, p.drag_kg_s, force0)
    after = np.maximum(time-p.empty_time_s, 0)
    tau = p.relaxation_time_s
    velocity = velocity_active*np.exp(-after/tau)
    position = position_active + velocity_active*tau*(-np.expm1(-after/tau))
    force = force0*ink_active(time, p)
    acceleration = (force-p.drag_kg_s*velocity)/p.mass_kg
    return velocity, acceleration, position


def simulate(p):
    """Exact linear-reservoir response, including coasting after ink depletion."""
    time = np.linspace(0, p.duration_s, p.points)
    # Include the depletion instant exactly, even on a coarse plotting grid.
    if 0 < p.empty_time_s <= p.duration_s:
        time = np.unique(np.append(time, [np.nextafter(p.empty_time_s, 0), p.empty_time_s]))
    fraction = ink_fraction(time, p)
    contrast = p.gamma_contrast_n_m*ink_active(time, p)
    drive = p.perimeter_m*contrast/4
    velocity, acceleration, position = finite_ink_response(time, p)
    drag_force = -p.drag_kg_s*velocity
    remaining = p.ink_volume_m3*fraction
    return {
        "time_s": time,
        "velocity_m_s": velocity,
        "acceleration_m_s2": acceleration,
        "position_m": position,
        "capillary_force_n": drive,
        "hydrodynamic_force_n": drag_force,
        "total_force_n": drive + drag_force,
        "contrast_n_m": contrast,
        "ink_remaining_m3": remaining,
        "ink_ejected_m3": p.ink_volume_m3-remaining,
        "shear_stress_pa": WATER_VISCOSITY_PA_S*velocity/p.shear_depth_m,
        "reynolds": WATER_DENSITY_KG_M3*np.abs(velocity)*(2*p.radius_m)/WATER_VISCOSITY_PA_S,
    }
