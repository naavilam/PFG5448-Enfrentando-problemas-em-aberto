"""Live controls: streamlit run app.py (requires Streamlit >= 1.37)."""
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from model import Params
from live_model import LiveSimulation

st.set_page_config(page_title='Ink Boat — controle ao vivo',layout='wide')
st.title('Ink Boat · Etapa 1 — controles ao vivo')
st.warning('Protótipo de filme fino sobre substrato rígido. Não representa explicitamente a interface tinta–água. Pétalas iniciais são perturbações impostas, não uma instabilidade validada.')
st.caption('Mude os controles durante a execução: o estado atual é preservado. Alterações na espessura inicial e resolução exigem reinicialização.')
with st.sidebar:
    st.header('Parâmetros ao vivo')
    mu=st.slider('Viscosidade μ (Pa·s)',0.005,0.5,0.08,0.005)
    rho=st.slider('Densidade ρ (kg/m³)',700,1500,1000,25)
    gamma=st.slider('Tensão superficial γ₀ (mN/m)',15.0,75.0,35.0,1.0)
    E=st.slider('Sensibilidade Marangoni EΓref (mN/m)',0.0,25.0,20.0,1.0)
    D=st.select_slider('Difusividade D (m²/s)',options=[1e-11,1e-10,1e-9,2e-9,1e-8,1e-7],value=2e-9)
    source=st.slider('Liberação relativa',0.0,2.0,0.6,0.1)
    anis=st.slider('Perturbação angular',0.0,0.4,0.05,0.01)
    petals=st.slider('Modo angular da fonte (n)',2,16,8)
    speed=st.slider('Avanço físico por atualização (s)',0.01,0.20,0.04,0.01)
    st.divider()
    h=st.slider('Espessura inicial (μm) — reiniciar',20,160,80,10)
    resolution=st.select_slider('Resolução — reiniciar',options=[48,64,96,128],value=64)
    a,b,c=st.columns(3)
    play=a.button('▶',help='Continuar')
    pause=b.button('⏸',help='Pausar')
    reset=c.button('↺',help='Reiniciar')

if 'live_engine' not in st.session_state or reset:
    st.session_state.live_engine=LiveSimulation(Params(n=resolution,thickness_um=h))
    st.session_state.live_running=False
if play: st.session_state.live_running=True
if pause: st.session_state.live_running=False

engine=st.session_state.live_engine
p=engine.p
p.mu_pa_s=mu; p.rho_kg_m3=rho; p.gamma_n_m=gamma*1e-3
p.elasticity_n_m=E*1e-3; p.diffusivity_m2_s=D
p.source_strength=source; p.anisotropy=anis; p.petals=petals
engine.refresh_source()

if (p.n!=resolution or p.thickness_um!=h):
    st.info('Resolução ou espessura inicial alterada: pressione ↺ para aplicar.')

@st.fragment(run_every='0.25s')
def live_view():
    eng=st.session_state.live_engine
    if st.session_state.live_running:
        advanced,steps=eng.advance(speed,max_steps=1000)
        if advanced<speed*0.99:
            st.warning('Passos numéricos muito pequenos: reduza resolução/avanço ou pause. A integração está limitada a 1000 passos por atualização.')
    st.metric('Tempo físico simulado',f'{eng.t:.3f} s')
    st.caption('Status: '+('▶ Executando' if st.session_state.live_running else '⏸ Pausado'))
    extent=[eng.x[0]*1e3,eng.x[-1]*1e3,eng.x[0]*1e3,eng.x[-1]*1e3]
    col1,col2=st.columns(2)
    with col1:
        fig,ax=plt.subplots(figsize=(5,4))
        im=ax.imshow(eng.G/1e-6,origin='lower',extent=extent,cmap='magma',vmin=0)
        fig.colorbar(im,ax=ax,label='Γ / Γref')
        ax.set(title='Espécie superficial',xlabel='x (mm)',ylabel='y (mm)')
        st.pyplot(fig);plt.close(fig)
    with col2:
        fig,ax=plt.subplots(figsize=(5,4))
        im=ax.imshow((eng.h-eng.h0)*1e6,origin='lower',extent=extent,cmap='coolwarm')
        fig.colorbar(im,ax=ax,label='Δh (μm)')
        ax.set(title='Deformação da película',xlabel='x (mm)',ylabel='y (mm)')
        st.pyplot(fig);plt.close(fig)

live_view()
