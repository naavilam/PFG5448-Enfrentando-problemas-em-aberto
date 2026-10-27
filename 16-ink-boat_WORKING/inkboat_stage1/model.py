"""Exploratory lubrication + insoluble surfactant model on a rigid substrate.
NOT a validated two-liquid model for ink on water.
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class Params:
    n: int = 96
    domain_mm: float = 30.0
    mu_pa_s: float = 0.08
    rho_kg_m3: float = 1000.0
    gamma_n_m: float = 0.035
    elasticity_n_m: float = 0.020
    diffusivity_m2_s: float = 2e-9
    thickness_um: float = 80.0
    source_strength: float = 0.6
    source_width_mm: float = 0.8
    anisotropy: float = 0.05
    petals: int = 8
    duration_s: float = 2.0
    frames: int = 24
    seed: int = 2


def run(p: Params):
    n=p.n
    L=p.domain_mm*1e-3
    dx=L/n
    x=(np.arange(n)-n/2)*dx
    xx,yy=np.meshgrid(x,x,indexing='xy')
    r=np.hypot(xx,yy)
    phi=np.arctan2(yy,xx)
    rng=np.random.default_rng(p.seed)
    h0=p.thickness_um*1e-6
    h=np.full((n,n),h0)
    # Insoluble surfactant surface concentration, normalized to a reference.
    G=np.zeros((n,n))
    perturb=1+p.anisotropy*np.cos(p.petals*phi)+0.015*rng.standard_normal((n,n))
    source=np.exp(-0.5*(r/(p.source_width_mm*1e-3))**2)*np.maximum(0.0,perturb)
    source/=source.sum()*dx*dx
    # A finite reservoir releases surface-active material over time.
    # Gamma_ref = 1e-6 kg/m2, source_strength is in Gamma_ref * m2/s units.
    Gref=1e-6
    total_mass_rate=p.source_strength*Gref*(np.pi*(1e-3)**2) # kg/s
    mu=p.mu_pa_s; rho=p.rho_kg_m3; gamma0=p.gamma_n_m
    E=p.elasticity_n_m/Gref # N m / kg (gamma=gamma0-E*Gamma)
    D=p.diffusivity_m2_s
    g=9.81
    def grad(a):
        return ((np.roll(a,-1,1)-np.roll(a,1,1))/(2*dx),
                (np.roll(a,-1,0)-np.roll(a,1,0))/(2*dx))
    def lap(a):
        return (np.roll(a,-1,0)+np.roll(a,1,0)+np.roll(a,-1,1)+np.roll(a,1,1)-4*a)/(dx*dx)
    def div(qx,qy):
        return (np.roll(qx,-1,1)-np.roll(qx,1,1)+np.roll(qy,-1,0)-np.roll(qy,1,0))/(2*dx)
    frames=[]
    # Conservative explicit stepping; adaptive timestep limited by fourth-order capillary term.
    t=0.0
    targets=np.linspace(0,p.duration_s,p.frames)
    for target in targets:
        while t < target-1e-12:
            maxh=max(float(np.max(h)),1e-9)
            dt_cap=0.025*3*mu*dx**4/(gamma0*maxh**3+1e-30)
            dt_mar=0.025*mu*dx**2/(max(p.elasticity_n_m,1e-9)*maxh+1e-30)
            dt_diff=0.15*dx*dx/(D+1e-20)
            dt=min(target-t,dt_cap,dt_mar,dt_diff,0.01)
            pressure=-gamma0*lap(h)+rho*g*(h-h0)
            px,py=grad(pressure)
            gx,gy=grad(G)
            # Marangoni stress = grad(gamma) = -E grad(G).
            tx,ty=-E*gx,-E*gy
            qx=-(h**3/(3*mu))*px+(h*h/(2*mu))*tx
            qy=-(h**3/(3*mu))*py+(h*h/(2*mu))*ty
            usx=-(h*h/(2*mu))*px+(h/mu)*tx
            usy=-(h*h/(2*mu))*py+(h/mu)*ty
            dh=-div(qx,qy)
            # Conservative advection-diffusion; central scheme may need upwinding
            # at large Pe. Clip below for prototype robustness only.
            dG=-div(G*usx,G*usy)+D*lap(G)+total_mass_rate*source
            h=np.maximum(h+dt*dh,h0*0.2)
            G=np.maximum(G+dt*dG,0.0)
            t+=dt
        frames.append((t,h.copy(),G.copy()))
    return x,frames
