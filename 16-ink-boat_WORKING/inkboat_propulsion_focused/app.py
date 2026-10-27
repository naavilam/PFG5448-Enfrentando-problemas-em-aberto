import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from model import Params,LiveSimulation

st.set_page_config(page_title='Ink Boat | Marangoni Propulsion',layout='wide')
st.markdown('<style>[data-testid="stSidebar"]{min-width:610px!important;width:610px!important;}@media(max-width:900px){[data-testid="stSidebar"]{min-width:350px!important;width:350px!important;}}</style>',unsafe_allow_html=True)
st.title('Ink Boat — Marangoni Propulsion')
st.markdown('**Propulsão do barco gerada por diferenças espaciais de tensão superficial**')
st.caption('Modelo: transporte superficial acoplado à resposta hidrodinâmica de Stokes.')
st.info('**O que é previsto:** força capilar e movimento translacional decorrentes da distribuição de tensão superficial calculada. **Limites:** Stokes (Re ≪ 1), superfície plana, fonte efetiva de tinta e arrasto aproximado; não resolve as três interfaces móveis nem as tensões sobre a superfície molhada.')

with st.sidebar:
    st.header('Modelo físico | Parâmetros')
    with st.container(border=True):
        st.markdown('**1 · Interfaces e contraste de tensão**')
        eq,ctrl=st.columns(2)
        with eq:
            st.latex(r'S=\gamma_{wa}-\gamma_{ia}-\gamma_{iw}')
            st.caption('S é diagnóstico: não fornece o gradiente de tensão.')
        with ctrl:
            wa=st.slider('γ água–ar (mN/m)',20.,80.,72.)
            ia=st.slider('γ tinta–ar (mN/m)',1.,70.,30.)
            iw=st.slider('γ tinta–água (mN/m)',1.,70.,25.)
        st.caption(f'S = {wa-ia-iw:+.2f} mN/m')
    with st.container(border=True):
        st.markdown('**2 · Como a tinta cria o gradiente de tensão**')
        eq,ctrl=st.columns(2)
        with eq:
            st.latex(r'\gamma=\gamma_{wa}-\beta\Gamma')
            st.latex(r'\partial_t\Gamma+\nabla_s\cdot(\Gamma\mathbf u_s)=D_s\nabla_s^2\Gamma+J-\Gamma/\tau')
            st.caption('Γ em kg/m²; J em kg/(m² s). β é uma hipótese constitutiva.')
        with ctrl:
            beta=st.slider('β (N m/kg)',10.,3000.,1000.,10.)
            ds=st.select_slider('Dₛ (m²/s)',options=[1e-10,1e-9,1e-8,1e-7,1e-6],value=1e-8,format_func=lambda x:f'{x:.0e}')
            release=st.slider('Liberação (ng/s)',0.,20.,1.,0.2)
            tau=st.slider('Tempo de remoção τ (s)',0.5,30.,8.,0.5)
            width=st.slider('Largura da fonte (mm)',0.5,5.,1.2,0.1)
            angle=st.slider('Direção da fonte (°)',0,359,180)
    with st.container(border=True):
        st.markdown('**3 · Escoamento induzido por Marangoni**')
        eq,ctrl=st.columns(2)
        with eq:
            st.latex(r'0=-\nabla p+\mu\nabla^2\mathbf u')
            st.latex(r'\nabla\cdot\mathbf u=0')
            st.latex(r'\widehat{\mathbf u}_s=\frac{(I-\hat k\hat k/2)\widehat{\nabla_s\gamma}}{\mu|k|}')
            st.caption('Solução espectral para água profunda, interface plana e domínio periódico.')
        with ctrl:
            mu=st.slider('μ água (mPa·s)',0.5,10.,1.,0.1)
            rho=st.slider('ρ água (kg/m³)',900.,1100.,998.)
            n=st.select_slider('Malha (reiniciar)',options=[48,64,96,128],value=96)
            domain=st.slider('Domínio (mm; reiniciar)',30.,100.,50.,5.)
    with st.container(border=True):
        st.markdown('**4 · Força propulsiva e velocidade**')
        eq,ctrl=st.columns(2)
        with eq:
            st.latex(r'\mathbf F=R\int_0^{2\pi}\gamma\mathbf e_r\,d\theta')
            st.latex(r'm\dot{\mathbf U}=\mathbf F-b\mathbf U')
            st.latex(r'b\approx 6\pi\mu R\,c_d')
            st.caption('b é aproximação de arrasto esférico, NÃO integração de tensões no barco.')
        with ctrl:
            radius=st.slider('Raio R (mm; reiniciar)',1.,8.,3.,0.5)
            mass=st.slider('Massa (g)',0.02,3.,0.2,0.02)
            drag=st.slider('Fator de arrasto c_d',0.3,15.,1.,0.1)
    dt=st.slider('Tempo simulado por atualização (s)',0.002,0.05,0.015,0.001)
    x,y,z=st.columns(3)
    start=x.button('▶ Iniciar',use_container_width=True)
    pause=y.button('⏸ Pausar',use_container_width=True)
    reset=z.button('↺ Reiniciar',use_container_width=True)

p=Params(n=n,domain_mm=domain,radius_mm=radius,mass_g=mass,mu_mpas=mu,rho=rho,
         gamma_wa_mnm=wa,gamma_ia_mnm=ia,gamma_iw_mnm=iw,
         beta_n_m_per_kg_m2=beta,ds=ds,release_ng_s=release,plume_width_mm=width,
         decay_s=tau,drag_factor=drag,source_angle_deg=angle)
if 'ink_coupled' not in st.session_state or reset:
    st.session_state.ink_coupled=LiveSimulation(p)
    st.session_state.ink_running=False
if start:st.session_state.ink_running=True
if pause:st.session_state.ink_running=False
s=st.session_state.ink_coupled
if (s.p.n,s.p.domain_mm,s.p.radius_mm)!=(n,domain,radius):
    st.info('Malha, domínio ou raio alterados: clique em Reiniciar para aplicar.')
for field in ('mass_g','mu_mpas','rho','gamma_wa_mnm','gamma_ia_mnm','gamma_iw_mnm',
              'beta_n_m_per_kg_m2','ds','release_ng_s','plume_width_mm','decay_s','drag_factor','source_angle_deg'):
    setattr(s.p,field,getattr(p,field))

@st.fragment(run_every='0.4s')
def live():
    s=st.session_state.ink_coupled
    if st.session_state.ink_running:
        s.step(dt)
    re,pe=s.diagnostics()
    # A interpretação principal é a força e o movimento do barco.
    speed=np.linalg.norm(s.vel)
    force=np.linalg.norm(s.force)
    a,b,c,d=st.columns(4)
    a.metric('Velocidade do barco',f'{speed*1e3:.3f} mm/s')
    b.metric('Força propulsiva',f'{force*1e6:.3f} µN')
    c.metric('Deslocamento',f'{np.linalg.norm(s.pos)*1e3:.3f} mm')
    d.metric('Tempo',f'{s.t:.2f} s')
    st.caption(f'Re = {re:.2g} (validade Stokes: Re ≪ 1) · Pe = {pe:.2g} · massa superficial = {np.sum(s.G)*s.dx**2*1e9:.3f} ng')
    left,right=st.columns([1.35,1.0],gap='large')
    with left:
        st.subheader('Resposta propulsiva')
        fig,ax=plt.subplots(figsize=(7,4.5))
        if s.history:
            h=np.asarray(s.history)
            ax.plot(h[:,0],np.hypot(h[:,3],h[:,4])*1e3,label='Velocidade |U| (mm/s)',lw=2)
            ax.set(xlabel='Tempo (s)',ylabel='Velocidade (mm/s)')
            ax.grid(alpha=.25)
        st.pyplot(fig,use_container_width=True);plt.close(fig)
        fig,ax=plt.subplots(figsize=(7,3.5))
        if s.history:
            h=np.asarray(s.history)
            ax.plot(h[:,0],np.hypot(h[:,5],h[:,6])*1e6,label='|F|',color='tab:orange',lw=2)
            ax.set(xlabel='Tempo (s)',ylabel='Força capilar (µN)')
            ax.grid(alpha=.25)
        st.pyplot(fig,use_container_width=True);plt.close(fig)
    with right:
        st.subheader('Mecanismo: tensão superficial')
        gamma=s.p.gamma_wa_mnm-s.p.beta_n_m_per_kg_m2*s.G*1e3
        fig,ax=plt.subplots(figsize=(5,4.2))
        extent=[-s.p.domain_mm/2,s.p.domain_mm/2]*2
        # Escala fixa: referência γ_wa e amplitude indicada pelo usuário.
        span=max(s.p.beta_n_m_per_kg_m2*1e-6*1e3,0.1)
        im=ax.imshow(gamma,extent=extent,origin='lower',cmap='viridis',
                     vmin=s.p.gamma_wa_mnm-10*span,vmax=s.p.gamma_wa_mnm)
        from matplotlib.patches import Circle
        ax.add_patch(Circle((s.pos[0]*1e3,s.pos[1]*1e3),s.p.radius_mm,
                            fill=False,color='white',lw=2))
        ax.set(xlabel='x (mm)',ylabel='y (mm)',title='γ água–ar efetiva (mN/m)',aspect='equal')
        fig.colorbar(im,ax=ax,label='mN/m')
        st.pyplot(fig,use_container_width=True);plt.close(fig)
        st.caption('O campo de γ é calculado a partir de Γ. A fonte angular de tinta ainda é uma entrada do modelo.')
        fig,ax=plt.subplots(figsize=(5,3.5))
        if s.history:
            h=np.asarray(s.history)
            ax.plot(h[:,1]*1e3,h[:,2]*1e3,lw=2)
            ax.plot(h[-1,1]*1e3,h[-1,2]*1e3,'o',ms=6)
        ax.set(xlabel='x (mm)',ylabel='y (mm)',title='Trajetória do barco',aspect='equal')
        ax.grid(alpha=.25)
        st.pyplot(fig,use_container_width=True);plt.close(fig)

live()
