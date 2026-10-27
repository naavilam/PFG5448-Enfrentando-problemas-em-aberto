"""Research prototype: insoluble surface species coupled to subphase Stokes flow.

Flat interface z=0, water occupies z<0, deep layer, Re<<1.
Fourier traction-to-velocity map for an incompressible 3D Stokes half-space:
 u_hat_s(k) = (I - kk/(2|k|^2)) tau_hat(k)/(mu |k|), k!=0.
This is NOT a 3-phase Navier-Stokes/VOF solver; boat drag is approximate.
"""
from dataclasses import dataclass
import numpy as np

@dataclass
class Params:
    n: int = 96
    domain_mm: float = 50.0
    radius_mm: float = 3.0
    mass_g: float = 0.20
    mu_mpas: float = 1.0
    rho: float = 998.0
    gamma_wa_mnm: float = 72.0
    gamma_ia_mnm: float = 30.0
    gamma_iw_mnm: float = 25.0
    beta_n_m_per_kg_m2: float = 1000.0
    ds: float = 1e-8
    release_ng_s: float = 1.0
    plume_width_mm: float = 1.2
    decay_s: float = 8.0
    drag_factor: float = 1.0
    source_angle_deg: float = 180.0

class LiveSimulation:
    def __init__(self, p=None):
        self.p = p or Params()
        p=self.p; self.dx=p.domain_mm*1e-3/p.n
        a=(np.arange(p.n)-p.n//2)*self.dx
        self.xx,self.yy=np.meshgrid(a,a)
        k=2*np.pi*np.fft.fftfreq(p.n,d=self.dx)
        self.kx,self.ky=np.meshgrid(k,k)
        self.k2=self.kx**2+self.ky**2
        self.k=np.sqrt(self.k2)
        self.G=np.zeros((p.n,p.n))
        self.pos=np.zeros(2);self.vel=np.zeros(2)
        self.t=0.;self.force=np.zeros(2);self.history=[]
        self.u=np.zeros_like(self.G);self.v=np.zeros_like(self.G)
        self.initial_mass=0.

    def grad(self,a):
        return ((np.roll(a,-1,axis=1)-np.roll(a,1,axis=1))/(2*self.dx),
                (np.roll(a,-1,axis=0)-np.roll(a,1,axis=0))/(2*self.dx))

    def source(self):
        p=self.p; R=p.radius_mm*1e-3
        ang=np.deg2rad(p.source_angle_deg)
        center=self.pos+R*np.array([np.cos(ang),np.sin(ang)])
        d2=(self.xx-center[0])**2+(self.yy-center[1])**2
        sig=max(p.plume_width_mm*1e-3,self.dx)
        s=np.exp(-d2/(2*sig**2))
        # Avoid depositing species inside the boat.
        r2=(self.xx-self.pos[0])**2+(self.yy-self.pos[1])**2
        s*=r2>=R**2
        return s/(np.sum(s)*self.dx**2+1e-40)

    def stokes_velocity(self):
        p=self.p; mu=p.mu_mpas*1e-3
        # Linear equation of state gamma=gamma_wa-beta*Gamma.
        gamma=p.gamma_wa_mnm*1e-3-p.beta_n_m_per_kg_m2*self.G
        tx,ty=self.grad(gamma)
        txh=np.fft.fft2(tx);tyh=np.fft.fft2(ty)
        k2=np.where(self.k2==0,1,self.k2)
        dot=(self.kx*txh+self.ky*tyh)/k2
        inv=1/np.where(self.k==0,np.inf,mu*self.k)
        uh=(txh-.5*self.kx*dot)*inv
        vh=(tyh-.5*self.ky*dot)*inv
        self.u=np.fft.ifft2(uh).real
        self.v=np.fft.ifft2(vh).real
        return gamma

    def contact_force(self,gamma):
        p=self.p;R=p.radius_mm*1e-3
        theta=np.linspace(0,2*np.pi,256,endpoint=False)
        xp=self.pos[0]+R*np.cos(theta)
        yp=self.pos[1]+R*np.sin(theta)
        # Bilinear periodic sampling of gamma at contact line.
        n=p.n
        gx=(xp/self.dx+n/2)%n; gy=(yp/self.dx+n/2)%n
        ix=np.floor(gx).astype(int);iy=np.floor(gy).astype(int)
        fx=gx-ix;fy=gy-iy
        g=((1-fx)*(1-fy)*gamma[iy,ix]+fx*(1-fy)*gamma[iy,(ix+1)%n]
           +(1-fx)*fy*gamma[(iy+1)%n,ix]+fx*fy*gamma[(iy+1)%n,(ix+1)%n])
        return R*(2*np.pi/len(theta))*np.array([np.dot(g,np.cos(theta)),np.dot(g,np.sin(theta))])

    def step(self,duration=0.02):
        p=self.p
        if p.n!=self.G.shape[0]:raise ValueError('Resolução alterada: reinicie.')
        left=float(duration)
        while left>1e-12:
            gamma=self.stokes_velocity()
            maxu=max(np.max(np.hypot(self.u,self.v)),np.linalg.norm(self.vel),1e-8)
            dt=min(left,0.15*self.dx/maxu,0.15*self.dx**2/max(p.ds,1e-30),0.025)
            # Conservative upwind finite-volume advection on a periodic domain.
            ux=.5*(self.u+np.roll(self.u,-1,axis=1))
            uy=.5*(self.v+np.roll(self.v,-1,axis=0))
            fluxx=ux*np.where(ux>0,self.G,np.roll(self.G,-1,axis=1))
            fluxy=uy*np.where(uy>0,self.G,np.roll(self.G,-1,axis=0))
            adv=-(fluxx-np.roll(fluxx,1,axis=1)+fluxy-np.roll(fluxy,1,axis=0))/self.dx
            lap=(np.roll(self.G,1,0)+np.roll(self.G,-1,0)+np.roll(self.G,1,1)+np.roll(self.G,-1,1)-4*self.G)/self.dx**2
            self.G=np.maximum(0,self.G+dt*(adv+p.ds*lap+p.release_ng_s*1e-12*self.source()-self.G/max(p.decay_s,1e-6)))
            gamma=p.gamma_wa_mnm*1e-3-p.beta_n_m_per_kg_m2*self.G
            self.force=self.contact_force(gamma)
            # Quasi-steady Stokes sphere drag is a placeholder for boat-fluid coupling.
            b=6*np.pi*p.mu_mpas*1e-3*p.radius_mm*1e-3*p.drag_factor
            mass=p.mass_g*1e-3
            vstar=self.force/b
            factor=np.exp(-b*dt/mass)
            old=self.vel.copy()
            self.vel=vstar+(old-vstar)*factor
            self.pos+=vstar*dt+(old-vstar)*(1-factor)*mass/b
            self.t+=dt;left-=dt
        self.history.append((self.t,*self.pos,*self.vel,*self.force,np.sum(self.G)*self.dx**2))
        if len(self.history)>3000:self.history=self.history[-3000:]

    def diagnostics(self):
        p=self.p; speed=np.linalg.norm(self.vel);L=2*p.radius_mm*1e-3
        return p.rho*speed*L/(p.mu_mpas*1e-3),speed*L/p.ds
