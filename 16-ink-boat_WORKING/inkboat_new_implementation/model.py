"""Equation 3.5 with uniform gamma and an explicit thin-disk shear closure.

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
    shear_depth_m: float = 0.5e-3
    duration_s: float = 5.0
    points: int = 501
    contour_points: int = 256
    external_force_n: float = 0.0  # optional validation load; not Marangoni propulsion

    def __post_init__(self):
        for name in ("mass_kg", "perimeter_m", "shear_depth_m", "duration_s"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not np.isfinite(self.gamma_n_m) or self.gamma_n_m < 0:
            raise ValueError("gamma_n_m must be finite and nonnegative")
        if not np.isfinite(self.external_force_n):
            raise ValueError("external_force_n must be finite")
        for name, minimum in (("points", 2), ("contour_points", 16)):
            value = getattr(self, name)
            if not isinstance(value, (int, np.integer)) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")

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

    Current model: one constant in space AND time. No ink transport is assumed.
    """
    return np.full_like(theta, p.gamma_n_m, dtype=float)


def capillary_force(p, gamma_values=None):
    """Closed-circle integral F = integral gamma*t dl (N).

    Removing a constant baseline is exact for a closed planar circle because
    integral t dl = 0. It prevents roundoff from manufacturing propulsion.
    Nonuniform arrays can be supplied for mathematical validation only;
    the application currently always calls gamma_on_contour (uniform).
    """
    theta, _, conormal = contact_line(p)
    gamma = gamma_on_contour(theta, 0.0, p) if gamma_values is None else np.asarray(gamma_values, dtype=float)
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


def simulate(p):
    """Evaluate the constant-gamma reduced model over the entire time window.

    Changing controls recomputes this window from U(0)=0, without blending
    trajectories from different parameter sets. No numerical time-step error.
    """
    time = np.linspace(0, p.duration_s, p.points)
    force_gamma = capillary_force(p)
    drive = float(force_gamma[0]) + p.external_force_n
    velocity, acceleration, position = constant_force_response(time, p.mass_kg, p.drag_kg_s, drive)
    drag_force = -p.drag_kg_s * velocity
    total = drive + drag_force
    return {
        "time_s": time,
        "velocity_m_s": velocity,
        "acceleration_m_s2": acceleration,
        "position_m": position,
        "capillary_force_n": np.full_like(time, force_gamma[0]),
        "hydrodynamic_force_n": drag_force,
        "external_force_n": np.full_like(time, p.external_force_n),
        "total_force_n": total,
        "shear_stress_pa": WATER_VISCOSITY_PA_S * velocity / p.shear_depth_m,
        "reynolds": WATER_DENSITY_KG_M3 * np.abs(velocity) * (2*p.radius_m) / WATER_VISCOSITY_PA_S,
    }
