import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from model import Params, LiveSimulation

st.set_page_config(page_title='Ink Boat — Propulsão Marangoni', layout='wide')
st.title('Ink Boat | Propulsão por tensão superficial')
st.info('Modelo reduzido em vista superior: calcula a força resultante na borda e a velocidade do barco. Não simula o espalhamento da tinta nem resolve Navier–Stokes multifásico.')
with st.sidebar:
    st.header('Interfaces (mN/m)')
    wa=st.slider('Água–ar γwa',20.,80.,72.,1.)
    ia=st.slider('Tinta–ar γia',10.,70.,30.,1.)
    iw=st.slider('Tinta–água γiw',1.,60.,25.,1.)
    st.header('Contraste efetivo na borda')
    dg=st.slider('Redução local Δγ (mN/m)',0.,40.,8.,0.5)
    release=st.slider('Intensidade da fonte (1/s)',0.,5.,0.8,0.1)
    angle=st.slider('Ângulo da fonte (graus)',0,359,180)
    width=st.slider('Largura angular (graus)',5,180,35)
    relax=st.slider('Tempo de relaxação (s)',0.1,20.,2.,0.1)
    diffusion=st.slider('Mistura angular (1/s)',0.,1.,0.15,0.01)
    st.header('Barco e resistência')
    radius=st.slider('Raio (mm)',1.,15.,4.,0.5)
    mass=st.slider('Massa (g)',0.05,5.,0.5,0.05)
    drag=st.slider('Arrasto linear b (g/s)',0.1,30.,4.,0.1)
    dt=st.slider('Avanço por atualização (s)',0.01,0.20,0.04,0.01)
    a,b,c=st.columns(3)
    play=a.button('▶'); pause=b.button('⏸'); reset=c.button('↺')

if 'engine' not in st.session_state or reset:
    st.session_state.engine=LiveSimulation()
    st.session_state.running=False
if play: st.session_state.running=True
if pause: st.session_state.running=False
p=st.session_state.engine.p
p.gamma_wa_mNm=wa;p.gamma_ia_mNm=ia;p.gamma_iw_mNm=iw
p.delta_gamma_mNm=min(dg,wa)
p.release_rate_s=release;p.release_angle_deg=angle;p.plume_width_deg=width
p.relaxation_s=relax;p.angular_diffusion_s=diffusion
p.radius_mm=radius;p.mass_g=mass;p.drag_kg_s=drag*1e-3

st.write(f'**Coeficiente de espalhamento termodinâmico:** S = γwa − γia − γiw = **{st.session_state.engine.spreading_coefficient_mNm:+.1f} mN/m**.')
st.caption('As três tensões são parâmetros físicos distintos. S indica a tendência termodinâmica de molhamento, mas NÃO determina sozinho o gradiente de tensão ao redor do barco. Δγ é um parâmetro efetivo independente, que requer calibração experimental.')

@st.fragment(run_every='0.3s')
def live():
    sim=st.session_state.engine
    if st.session_state.running: sim.step(dt)
    fx,fy=sim.net_force()
    c1,c2,c3,c4=st.columns(4)
    c1.metric('Tempo',f'{sim.t:.2f} s')
    c2.metric('Força longitudinal',f'{fx*1e6:.2f} μN')
    c3.metric('Velocidade',f'{sim.v*1e3:.2f} mm/s')
    c4.metric('Deslocamento',f'{sim.x*1e3:.2f} mm')
    left,right=st.columns(2)
    with left:
        fig,ax=plt.subplots(figsize=(5,4))
        theta=sim.theta
        gamma=sim.gamma_profile()*1e3
        ax.plot(np.rad2deg(theta),gamma)
        ax.set(xlabel='Ângulo na borda (graus)',ylabel='γ efetiva água–ar (mN/m)',title='Tensão superficial ao redor do barco',xlim=(0,360),ylim=(max(0,wa-dg-2),wa+2))
        ax.grid(alpha=0.3)
        st.pyplot(fig);plt.close(fig)
    with right:
        fig,ax=plt.subplots(figsize=(5,4))
        if sim.history:
            hist=np.asarray(sim.history)
            ax.plot(hist[:,0],hist[:,2]*1e3,label='v (mm/s)')
            ax.set(xlabel='Tempo (s)',ylabel='Velocidade (mm/s)',title='Resposta dinâmica do barco')
            ax.grid(alpha=0.3)
        st.pyplot(fig);plt.close(fig)
    st.caption('F = R ∫ γ(θ) (cos θ, sin θ) dθ; m dv/dt = Fx − b v. O modelo assume movimento apenas ao longo de x, arrasto linear e uma distribuição prescrita/relaxante de tensão superficial.')
live()
