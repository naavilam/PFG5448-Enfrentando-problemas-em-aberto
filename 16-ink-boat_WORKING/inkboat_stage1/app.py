"""Run: streamlit run app.py"""
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from model import Params,run

st.set_page_config(page_title='Ink Boat — Etapa 1',layout='wide')
st.title('Ink Boat · Etapa 1 — espalhamento com Marangoni')
st.warning('Protótipo exploratório: filme fino sobre substrato rígido + espécie superficial insolúvel. NÃO é uma simulação validada de tinta imiscível sobre água. Estrelas/pétalas não são garantidas pela física deste modelo.')
st.markdown('A tensão superficial varia com a concentração superficial: **γ = γ₀ − EΓ**. O gradiente de γ produz uma tensão de Marangoni, que transporta o material e altera a espessura do filme.')
with st.sidebar:
    st.header('Parâmetros')
    mu=st.slider('Viscosidade μ (Pa·s)',0.005,0.5,0.08,0.005)
    rho=st.slider('Densidade ρ (kg/m³)',700,1500,1000,25)
    gamma=st.slider('Tensão superficial de referência γ₀ (mN/m)',15.0,75.0,35.0,1.0)
    E=st.slider('Sensibilidade Marangoni EΓref (mN/m)',0.0,25.0,20.0,1.0)
    D=st.select_slider('Difusividade superficial D (m²/s)',options=[1e-11,1e-10,1e-9,2e-9,1e-8,1e-7],value=2e-9)
    h=st.slider('Espessura inicial do filme (μm)',20,160,80,10)
    source=st.slider('Liberação relativa de espécie ativa',0.0,2.0,0.6,0.1)
    anis=st.slider('Perturbação angular inicial (amplitude)',0.0,0.4,0.05,0.01)
    petals=st.slider('Modo angular inicial (n)',2,16,8)
    duration=st.slider('Tempo simulado (s)',0.2,5.0,2.0,0.2)
    resolution=st.select_slider('Resolução',options=[48,64,96,128],value=64)
    start=st.button('Simular',type='primary')
if start:
    p=Params(n=resolution,mu_pa_s=mu,rho_kg_m3=rho,gamma_n_m=gamma*1e-3,
             elasticity_n_m=E*1e-3,diffusivity_m2_s=D,thickness_um=h,
             source_strength=source,anisotropy=anis,petals=petals,duration_s=duration)
    with st.spinner('Resolvendo as equações...'):
        x,frames=run(p)
    selected=st.slider('Instante',0,len(frames)-1,len(frames)-1)
    t,H,G=frames[selected]
    c1,c2=st.columns(2)
    extent=[x[0]*1e3,x[-1]*1e3,x[0]*1e3,x[-1]*1e3]
    with c1:
        fig,ax=plt.subplots(figsize=(5,4))
        im=ax.imshow(G/1e-6,origin='lower',extent=extent,cmap='magma')
        fig.colorbar(im,ax=ax,label='Γ / Γref')
        ax.set(title=f'Espécie superficial · t={t:.2f} s',xlabel='x (mm)',ylabel='y (mm)')
        st.pyplot(fig);plt.close(fig)
    with c2:
        fig,ax=plt.subplots(figsize=(5,4))
        im=ax.imshow((H-p.thickness_um*1e-6)*1e6,origin='lower',extent=extent,cmap='coolwarm')
        fig.colorbar(im,ax=ax,label='Δh (μm)')
        ax.set(title='Deformação da película',xlabel='x (mm)',ylabel='y (mm)')
        st.pyplot(fig);plt.close(fig)
    st.caption('Contornos e estruturas são resultados do modelo reduzido; a perturbação angular inicial é imposta pelo controle. O campo de Γ não é a fronteira de uma fase líquida distinta.')
else:
    st.info('Ajuste os parâmetros e clique em **Simular**.')
