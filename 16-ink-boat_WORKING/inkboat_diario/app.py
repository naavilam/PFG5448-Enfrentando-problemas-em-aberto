import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from model import LiveSimulation

st.set_page_config(page_title='Ink Boat | Modelo do diário', layout='wide')
st.markdown('''<style>[data-testid="stSidebar"] {min-width: 600px !important; width: 600px !important;} @media (max-width: 900px) {[data-testid="stSidebar"] {min-width: 340px !important; width: 340px !important;}}</style>''', unsafe_allow_html=True)
st.title('Ink Boat — Propulsão de Marangoni')
st.warning('Modelo simplificado: incorpora os balanços de força e torque, mas **não resolve Navier–Stokes 3D, transporte de tinta no volume nem tensões hidrodinâmicas na superfície molhada**. O arrasto e a distribuição de tensão são efetivos. As três tensões não bastam para prever a força.')

with st.sidebar:
    st.header('Equações, constantes e parâmetros')
    with st.container(border=True):
        st.markdown('**1 · Interfaces**')
        eq, ct=st.columns(2)
        with eq:
            st.latex(r'S=\gamma_{wa}-\gamma_{ia}-\gamma_{iw}')
            st.caption('S descreve molhamento; não é a força propulsiva.')
        with ct:
            wa=st.slider('Água–ar (mN/m)',20.,80.,72.)
            ia=st.slider('Tinta–ar (mN/m)',1.,70.,30.)
            iw=st.slider('Tinta–água (mN/m)',1.,70.,25.)
        st.caption(f'S = {wa-ia-iw:+.2f} mN/m')
    with st.container(border=True):
        st.markdown('**2 · Fonte e transporte angular efetivo**')
        eq,ct=st.columns(2)
        with eq:
            st.latex(r'q(\theta)=\exp[-d_{2\pi}^2/(2\sigma^2)]')
            st.latex(r'\partial_t c=kq(1-c)-c/\tau+D_\theta\partial_{\theta\theta}c')
            st.caption('c(θ,0)=0; N=256. Transporte apenas na borda, NÃO o transporte 3D do artigo.')
        with ct:
            angle=st.slider('Direção da fonte (°)',0,359,180)
            width=st.slider('Largura (°)',5,180,35)
            rate=st.slider('Taxa k (s⁻¹)',0.,5.,0.8,0.1)
            tau=st.slider('Relaxação τ (s)',0.1,20.,2.,0.1)
            diff=st.slider('Mistura angular Dθ (s⁻¹)',0.,1.,0.15,0.01)
    with st.container(border=True):
        st.markdown('**3 · Tensão e força de linha de contato**')
        eq,ct=st.columns(2)
        with eq:
            st.latex(r'\gamma=\gamma_{wa}-\Delta\gamma\,c')
            st.latex(r'\mathbf F_{st}=R\int_0^{2\pi}\gamma\,\mathbf e_r\,d\theta')
            st.caption('A contribuição de tensão viscosa na parte submersa é aproximada pelo arrasto.')
        with ct:
            dg=st.slider('Contraste Δγ (mN/m)',0.,40.,8.,0.5)
            radius=st.slider('Raio R (mm)',1.,15.,4.,0.5)
    with st.container(border=True):
        st.markdown('**4 · Translação (x, y)**')
        eq,ct=st.columns(2)
        with eq:
            st.latex(r'm\dot{\mathbf U}=\mathbf F_{st}-b\mathbf U')
            st.latex(r'\dot{\mathbf X}=\mathbf U')
            st.caption('X(0)=U(0)=0. Arrasto linear substitui integral de tensões sobre a superfície molhada.')
        with ct:
            mass=st.slider('Massa m (g)',0.05,5.,0.5,0.05)
            drag=st.slider('Arrasto b (g/s)',0.1,30.,4.,0.1)
    with st.container(border=True):
        st.markdown('**5 · Rotação e torque**')
        eq,ct=st.columns(2)
        with eq:
            st.latex(r'I_{zz}\dot\omega=\tau_z-b_r\omega')
            st.latex(r'\dot\phi=\omega')
            st.latex(r'\tau_z=(\mathbf r_0\times\mathbf F_{st})_z')
            st.caption('Izz=κmR². Para borda circular com tração radial, torque é zero. Braço excêntrico é hipótese adicional.')
        with ct:
            inertia=st.slider('Fator de inércia κ',0.1,1.,0.5,0.05)
            offset=st.slider('Braço lateral efetivo (mm)',-4.,4.,0.,0.1)
            rot_drag=st.slider('Arrasto rotacional (10⁻⁹ N m s)',0.1,20.,1.,0.1)
    with st.container(border=True):
        st.markdown('**6 · Regime e integração**')
        eq,ct=st.columns(2)
        with eq:
            st.latex(r'Re=\rho U(2R)/\mu')
            st.latex(r'Pe=U(2R)/D_s')
            st.caption('Pe é apenas diagnóstico: D_s é uma hipótese para espécie interfacial, não difusão de tinta dissolvida.')
        with ct:
            rho=st.slider('Densidade água (kg/m³)',900.,1100.,1000.)
            mu=st.slider('Viscosidade água (mPa·s)',0.5,5.,1.,0.1)
            ds=st.select_slider('Difusividade superficial (m²/s)',options=[1e-11,1e-10,1e-9,1e-8,1e-7],value=1e-9,format_func=lambda v:f'{v:.0e}')
            dt=st.slider('Tempo por atualização (s)',0.01,0.20,0.04,0.01)
    a,b,c=st.columns(3)
    start=a.button('▶ Iniciar',use_container_width=True)
    stop=b.button('⏸ Pausar',use_container_width=True)
    reset=c.button('↺ Reiniciar',use_container_width=True)

if 'sim' not in st.session_state or reset:
    st.session_state.sim=LiveSimulation()
    st.session_state.running=False
if start: st.session_state.running=True
if stop: st.session_state.running=False
sim=st.session_state.sim
p=sim.p
p.gamma_wa_mNm,p.gamma_ia_mNm,p.gamma_iw_mNm=wa,ia,iw
p.delta_gamma_mNm=min(dg,wa)
p.release_angle_deg,p.plume_width_deg=angle,width
p.release_rate_s,p.relaxation_s,p.angular_diffusion_s=rate,tau,diff
p.radius_mm,p.mass_g,p.drag_g_s=radius,mass,drag
p.inertia_factor,p.source_offset_mm,p.rotational_drag_nms=inertia,offset,rot_drag*1e-9
p.rho_water,p.mu_water_mPas,p.surface_diffusivity_m2s=rho,mu,ds

@st.fragment(run_every='0.3s')
def display():
    s=st.session_state.sim
    if st.session_state.running: s.step(dt)
    re,pe=s.dimensionless()
    a,b,c,d=st.columns(4)
    a.metric('Tempo',f'{s.t:.2f} s')
    b.metric('Velocidade',f'{np.linalg.norm(s.vel)*1e3:.3f} mm/s')
    c.metric('Força',f'{np.linalg.norm(s.force)*1e6:.3f} μN')
    d.metric('Deslocamento',f'{np.linalg.norm(s.pos)*1e3:.3f} mm')
    st.caption(f'Re = {re:.3g} · Pe = {pe:.3g} · ângulo = {np.rad2deg(s.heading):.2f}° · ω = {s.omega:.3g} rad/s · torque = {s.torque*1e9:.3g} nN·m')
    left,right=st.columns(2)
    with left:
        fig,ax=plt.subplots(figsize=(6,4))
        ax.plot(np.rad2deg(s.theta),s.gamma_profile()*1e3)
        ax.set(xlim=(0,360),ylim=(max(0,wa-min(dg,wa)-2),wa+2),xlabel='Ângulo no barco (°)',ylabel='γ água–ar efetiva (mN/m)',title='Tensão na linha de contato')
        ax.grid(alpha=.3); st.pyplot(fig);plt.close(fig)
    with right:
        fig,ax=plt.subplots(figsize=(6,4))
        if s.history:
            h=np.array(s.history)
            ax.plot(h[:,0],h[:,5]*1e3,label='|U|')
            ax.legend()
        ax.set(xlabel='Tempo (s)',ylabel='Velocidade (mm/s)',title='Resposta translacional')
        ax.grid(alpha=.3);st.pyplot(fig);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,3))
    if s.history:
        h=np.array(s.history)
        ax.plot(h[:,1]*1e3,h[:,2]*1e3,lw=2)
        ax.scatter([h[-1,1]*1e3],[h[-1,2]*1e3],s=45)
    ax.set(xlabel='x (mm)',ylabel='y (mm)',title='Trajetória (vista superior)')
    ax.axis('equal');ax.grid(alpha=.3);st.pyplot(fig);plt.close(fig)
    st.caption('A formulação completa do diário exige Navier–Stokes, condição de cisalhamento Marangoni e transporte de espécie acoplado. Esta versão é um modelo de corpo rígido com fechamento fenomenológico.')

display()
