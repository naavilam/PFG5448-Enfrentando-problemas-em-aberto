"""Reduced 2-D (top-view) Marangoni propulsion model.

Prescribes the surface-tension contrast along a circular boat perimeter;
no hydrodynamic flow or three-phase contact-line problem is solved.
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class Params:
    radius_mm: float = 4.0
    mass_g: float = 0.5
    drag_kg_s: float = 0.004
    gamma_wa_mNm: float = 72.0
    gamma_ia_mNm: float = 30.0
    gamma_iw_mNm: float = 25.0
    delta_gamma_mNm: float = 8.0
    release_rate_s: float = 0.8
    relaxation_s: float = 2.0
    angular_diffusion_s: float = 0.15
    plume_width_deg: float = 35.0
    release_angle_deg: float = 180.0
    n: int = 256

class LiveSimulation:
    def __init__(self, p=None):
        self.p = p or Params()
        self.theta = 2*np.pi*np.arange(self.p.n)/self.p.n
        self.dtheta = 2*np.pi/self.p.n
        self.c = np.zeros(self.p.n)
        self.x = 0.0
        self.v = 0.0
        self.t = 0.0
        self.force = 0.0
        self.history = []

    def gamma_profile(self):
        p = self.p
        # Effective water-air surface tension adjacent to the hull.
        return (p.gamma_wa_mNm - p.delta_gamma_mNm*self.c)*1e-3

    def net_force(self):
        # F = R ∮ gamma(theta) n(theta) dtheta; uniform gamma integrates to zero.
        g = self.gamma_profile()
        R = self.p.radius_mm*1e-3
        fx = R*self.dtheta*np.sum(g*np.cos(self.theta))
        fy = R*self.dtheta*np.sum(g*np.sin(self.theta))
        return fx, fy

    def step(self, duration=0.03):
        p = self.p
        width = max(np.deg2rad(p.plume_width_deg), 0.01)
        center = np.deg2rad(p.release_angle_deg)
        angular_dist = np.angle(np.exp(1j*(self.theta-center)))
        source = np.exp(-0.5*(angular_dist/width)**2)
        # Bounded concentration, 0 <= c <= 1.  Units are phenomenological.
        remaining = duration
        while remaining > 1e-12:
            dt = min(remaining, 0.01, 0.20*self.dtheta**2/max(p.angular_diffusion_s,1e-12))
            lap = (np.roll(self.c,-1)-2*self.c+np.roll(self.c,1))/self.dtheta**2
            dc = p.release_rate_s*source*(1-self.c) - self.c/max(p.relaxation_s,1e-6) + p.angular_diffusion_s*lap
            self.c = np.clip(self.c + dt*dc, 0, 1)
            fx, _ = self.net_force()
            self.force = fx
            m = p.mass_g*1e-3
            b = p.drag_kg_s
            # Exact integration for constant force over the small time step.
            v_eq = fx/b
            decay = np.exp(-b*dt/m)
            v_old = self.v
            self.v = v_eq + (v_old-v_eq)*decay
            self.x += v_eq*dt + (v_old-v_eq)*(1-decay)*m/b
            self.t += dt
            remaining -= dt
        self.history.append((self.t, self.x, self.v, self.force))
        if len(self.history)>3000:
            self.history=self.history[-3000:]

    @property
    def spreading_coefficient_mNm(self):
        p=self.p
        return p.gamma_wa_mNm-p.gamma_ia_mNm-p.gamma_iw_mNm
