"""Incremental integration of the original exploratory thin-film model.
This is NOT a validated ink-water free-surface model.
"""
import numpy as np
from model import Params

class LiveSimulation:
    def __init__(self, p: Params):
        self.p=p
        n=p.n; self.dx=p.domain_mm*1e-3/n
        self.x=(np.arange(n)-n/2)*self.dx
        xx,yy=np.meshgrid(self.x,self.x,indexing='xy')
        self.r=np.hypot(xx,yy); self.phi=np.arctan2(yy,xx)
        self.noise=np.random.default_rng(p.seed).standard_normal((n,n))
        self.h0=p.thickness_um*1e-6
        self.h=np.full((n,n),self.h0)
        self.G=np.zeros((n,n)); self.t=0.0
        self.refresh_source()
    def refresh_source(self):
        p=self.p
        perturb=1+p.anisotropy*np.cos(p.petals*self.phi)+0.015*self.noise
        source=np.exp(-0.5*(self.r/(p.source_width_mm*1e-3))**2)*np.maximum(0,perturb)
        self.source=source/(source.sum()*self.dx**2)
    def advance(self, physical_seconds=0.03, max_steps=1000):
        p=self.p; dx=self.dx; mu=p.mu_pa_s; rho=p.rho_kg_m3
        gamma0=p.gamma_n_m; E=p.elasticity_n_m/1e-6; D=p.diffusivity_m2_s
        h=self.h; G=self.G; h0=self.h0
        rate=p.source_strength*1e-6*np.pi*1e-6
        def grad(a):
            return ((np.roll(a,-1,1)-np.roll(a,1,1))/(2*dx),
                    (np.roll(a,-1,0)-np.roll(a,1,0))/(2*dx))
        def lap(a):
            return (np.roll(a,-1,0)+np.roll(a,1,0)+np.roll(a,-1,1)+np.roll(a,1,1)-4*a)/dx**2
        def div(x,y):
            return (np.roll(x,-1,1)-np.roll(x,1,1)+np.roll(y,-1,0)-np.roll(y,1,0))/(2*dx)
        remaining=physical_seconds; steps=0
        while remaining>1e-12 and steps<max_steps:
            maxh=max(float(np.max(h)),1e-9)
            dt_cap=0.025*3*mu*dx**4/(gamma0*maxh**3+1e-30)
            dt_mar=0.025*mu*dx**2/(max(p.elasticity_n_m,1e-9)*maxh+1e-30)
            dt_diff=0.15*dx**2/(D+1e-20)
            dt=min(remaining,dt_cap,dt_mar,dt_diff,0.01)
            if dt<1e-10: break
            pressure=-gamma0*lap(h)+rho*9.81*(h-h0)
            px,py=grad(pressure); gx,gy=grad(G)
            tx,ty=-E*gx,-E*gy
            qx=-(h**3/(3*mu))*px+h*h/(2*mu)*tx
            qy=-(h**3/(3*mu))*py+h*h/(2*mu)*ty
            usx=-(h*h/(2*mu))*px+h/mu*tx
            usy=-(h*h/(2*mu))*py+h/mu*ty
            dh=-div(qx,qy)
            dG=-div(G*usx,G*usy)+D*lap(G)+rate*self.source
            h=np.maximum(h+dt*dh,h0*0.2)
            G=np.maximum(G+dt*dG,0)
            self.t+=dt; remaining-=dt; steps+=1
        self.h=h; self.G=G
        return physical_seconds-remaining,steps
