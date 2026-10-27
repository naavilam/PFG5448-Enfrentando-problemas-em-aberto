"""Ink Boat: reduced 2-D rigid-body model inspired by Kang's force/torque balances.
NOT a Navier–Stokes or three-phase flow solver. The interfacial concentration is a
phenomenological field on the contact perimeter, not dissolved bulk concentration.
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class Params:
    radius_mm: float = 4.0
    mass_g: float = 0.5
    inertia_factor: float = 0.5  # Izz = factor*m*R^2
    drag_g_s: float = 4.0
    rotational_drag_nms: float = 1e-9
    gamma_wa_mNm: float = 72.0
    gamma_ia_mNm: float = 30.0
    gamma_iw_mNm: float = 25.0
    delta_gamma_mNm: float = 8.0
    release_rate_s: float = 0.8
    relaxation_s: float = 2.0
    angular_diffusion_s: float = 0.15
    plume_width_deg: float = 35.0
    release_angle_deg: float = 180.0
    source_offset_mm: float = 0.0
    rho_water: float = 1000.0
    mu_water_mPas: float = 1.0
    surface_diffusivity_m2s: float = 1e-9
    n: int = 256

class LiveSimulation:
    def __init__(self, p=None):
        self.p = p or Params()
        self.theta = 2*np.pi*np.arange(self.p.n)/self.p.n
        self.dtheta = 2*np.pi/self.p.n
        self.c = np.zeros(self.p.n)
        self.pos = np.zeros(2)
        self.vel = np.zeros(2)
        self.heading = 0.0
        self.omega = 0.0
        self.t = 0.0
        self.force = np.zeros(2)
        self.torque = 0.0
        self.history = []

    def gamma_profile(self):
        p = self.p
        return 1e-3*(p.gamma_wa_mNm - p.delta_gamma_mNm*self.c)

    def source_profile(self):
        p = self.p
        center = np.deg2rad(p.release_angle_deg)
        d = np.angle(np.exp(1j*(self.theta-center)))
        return np.exp(-0.5*(d/max(np.deg2rad(p.plume_width_deg), 1e-4))**2)

    def capillary_force_body(self):
        # Uniform gamma cancels exactly for the evenly sampled circle.
        R = self.p.radius_mm*1e-3
        unit = np.stack((np.cos(self.theta), np.sin(self.theta)), axis=1)
        return R*self.dtheta*np.sum(self.gamma_profile()[:,None]*unit, axis=0)

    def step(self, duration=0.04):
        p = self.p
        R = p.radius_mm*1e-3
        m = p.mass_g*1e-3
        I = p.inertia_factor*m*R*R
        b = p.drag_g_s*1e-3
        br = p.rotational_drag_nms
        source = self.source_profile()
        remaining = duration
        while remaining > 1e-12:
            dt = min(remaining, 0.01, 0.20*self.dtheta**2/max(p.angular_diffusion_s,1e-12))
            lap = (np.roll(self.c,-1)-2*self.c+np.roll(self.c,1))/self.dtheta**2
            # Angular advection in boat coordinates is NOT added: the field is
            # anchored to the rotating boat, and water flow is not solved.
            dc = p.release_rate_s*source*(1-self.c)-self.c/max(p.relaxation_s,1e-9)+p.angular_diffusion_s*lap
            self.c = np.clip(self.c+dt*dc, 0.0, 1.0)
            fb = self.capillary_force_body()
            ch, sh = np.cos(self.heading), np.sin(self.heading)
            rot = np.array([[ch,-sh],[sh,ch]])
            self.force = rot@fb
            # For a circular contact line, radial capillary tractions give ZERO
            # torque. An eccentric force application is an explicit, optional
            # phenomenological approximation, not predicted by circular geometry.
            offset_body = np.array([0.0, p.source_offset_mm*1e-3])
            arm = rot@offset_body
            self.torque = arm[0]*self.force[1]-arm[1]*self.force[0]
            v0 = self.vel.copy()
            decay = np.exp(-b*dt/m)
            veq = self.force/b
            self.vel = veq+(v0-veq)*decay
            self.pos += veq*dt+(v0-veq)*(1-decay)*m/b
            w0 = self.omega
            decay_r = np.exp(-br*dt/I)
            weq = self.torque/br
            self.omega = weq+(w0-weq)*decay_r
            self.heading += weq*dt+(w0-weq)*(1-decay_r)*I/br
            self.t += dt
            remaining -= dt
        self.history.append((self.t,*self.pos,*self.vel,np.linalg.norm(self.vel),*self.force,self.heading,self.omega,self.torque))
        if len(self.history)>4000:
            self.history=self.history[-4000:]

    @property
    def spreading_coefficient_mNm(self):
        p=self.p
        return p.gamma_wa_mNm-p.gamma_ia_mNm-p.gamma_iw_mNm

    def dimensionless(self):
        p=self.p
        U=np.linalg.norm(self.vel)
        L=2*p.radius_mm*1e-3
        re=p.rho_water*U*L/(p.mu_water_mPas*1e-3)
        pe=U*L/max(p.surface_diffusivity_m2s,1e-30)
        return re,pe
