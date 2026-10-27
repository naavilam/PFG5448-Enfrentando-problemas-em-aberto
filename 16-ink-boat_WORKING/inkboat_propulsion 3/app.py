"""Ink Boat: controls and governing equations paired in the sidebar."""
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from model import LiveSimulation

st.set_page_config(page_title='Ink Boat — Propulsão Marangoni', layout='wide')
st.markdown('''<style>
[data-testid="stSidebar"] {min-width: 590px !important; width: 590px !important;}
[data-testid="stSidebar"] [data-testid="stVerticalBlock"] {gap: .45rem;}
[data-testid="stSidebar"] .stSlider {padding-top: 0;}
@media (max-width: 1000px) { [data-testid="stSidebar"] {min-width: 340px !important; width: 340px !important;} }
</style>''', unsafe_allow_html=True)
st.title('Ink Boat | Propulsão por tensão superficial')
st.info('Modelo reduzido de propulsão: a distribuição de tensão é prescrita por uma fonte angular; não resolvemos o escoamento multifásico nem o espalhamento da tinta.')

with st.sidebar:
    st.title('Modelo e parâmetros')
    st.caption('Em cada bloco: equações e constantes à esquerda; parâmetros ajustáveis à direita.')

    with st.container(border=True):
        st.markdown('**1 · Três interfaces e molhamento**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'S=\gamma_{wa}-\gamma_{ia}-\gamma_{iw}')
            st.caption('S: coeficiente de espalhamento (mN/m). Não determina, sozinho, a propulsão.')
        with controls:
            wa = st.slider('Água–ar γwa (mN/m)', 20., 80., 72., 1.)
            ia = st.slider('Tinta–ar γia (mN/m)', 10., 70., 30., 1.)
            iw = st.slider('Tinta–água γiw (mN/m)', 1., 60., 25., 1.)
        st.caption(f'**S = {wa - ia - iw:+.1f} mN/m**')

    with st.container(border=True):
        st.markdown('**2 · Fonte angular**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'q(\theta)=e^{-d_{2\pi}(\theta,\theta_0)^2/(2\sigma_\theta^2)}')
            st.caption('d₂π: menor distância angular periódica. Constantes: 2π rad = 360°; amplitude máxima de q = 1.')
        with controls:
            angle = st.slider('Direção θ₀ (°)', 0, 359, 180)
            width = st.slider('Largura σθ (°)', 5, 180, 35)

    with st.container(border=True):
        st.markdown('**3 · Evolução da espécie superficial**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'\partial_t c=kq(1-c)-c/\tau+D_\theta\partial_{\theta\theta}c')
            st.caption('c: fração adimensional, 0 ≤ c ≤ 1; c(θ,0) = 0. Dθ: mistura angular fenomenológica, não difusividade molecular.')
        with controls:
            release = st.slider('Fonte k (s⁻¹)', 0., 5., 0.8, 0.1)
            relax = st.slider('Relaxação τ (s)', 0.1, 20., 2., 0.1)
            diffusion = st.slider('Mistura Dθ (s⁻¹)', 0., 1., 0.15, 0.01)

    with st.container(border=True):
        st.markdown('**4 · Tensão superficial efetiva**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'\gamma(\theta,t)=\gamma_{wa}-\Delta\gamma\,c(\theta,t)')
            st.caption('γ em N/m no cálculo; 1 mN/m = 10⁻³ N/m. Δγ é independente de S e requer calibração.')
        with controls:
            dg = st.slider('Contraste Δγ (mN/m)', 0., 40., 8., 0.5)
            if dg > wa:
                st.warning('Δγ foi limitado a γwa para evitar γ < 0.')

    with st.container(border=True):
        st.markdown('**5 · Força capilar resultante**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'F_x=R\int_0^{2\pi}\gamma(\theta,t)\cos\theta\,d\theta')
            st.latex(r'F_y=R\int_0^{2\pi}\gamma(\theta,t)\sin\theta\,d\theta')
            st.caption('Integração na borda circular; 1 mm = 10⁻³ m. O movimento usa apenas Fx.')
        with controls:
            radius = st.slider('Raio R (mm)', 1., 15., 4., 0.5)

    with st.container(border=True):
        st.markdown('**6 · Movimento e arrasto**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'm\dot v=F_x-bv')
            st.latex(r'\dot x=v')
            st.caption('Condições iniciais: x(0)=0, v(0)=0. Conversões: 1 g = 10⁻³ kg; 1 g/s = 10⁻³ kg/s.')
        with controls:
            mass = st.slider('Massa m (g)', 0.05, 5., 0.5, 0.05)
            drag = st.slider('Arrasto b (g/s)', 0.1, 30., 4., 0.1)

    with st.container(border=True):
        st.markdown('**7 · Integração numérica**')
        eq, controls = st.columns([1, 1.05], gap='medium')
        with eq:
            st.latex(r'\Delta\theta=2\pi/N')
            st.caption('N = 256 pontos angulares; passo interno ≤ 0,01 s, reduzido por estabilidade da mistura.')
        with controls:
            dt = st.slider('Avanço por atualização (s)', 0.01, 0.20, 0.04, 0.01)
        play_col, pause_col, reset_col = st.columns(3)
        play = play_col.button('▶ Iniciar', use_container_width=True)
        pause = pause_col.button('⏸ Pausar', use_container_width=True)
        reset = reset_col.button('↺ Reiniciar', use_container_width=True)

    st.caption('Modelo efetivo: não calcula hidrodinâmica de três fases nem a evolução da linha de contato.')

if 'engine' not in st.session_state or reset:
    st.session_state.engine = LiveSimulation()
    st.session_state.running = False
if play:
    st.session_state.running = True
if pause:
    st.session_state.running = False

p = st.session_state.engine.p
p.gamma_wa_mNm = wa
p.gamma_ia_mNm = ia
p.gamma_iw_mNm = iw
p.delta_gamma_mNm = min(dg, wa)
p.release_rate_s = release
p.release_angle_deg = angle
p.plume_width_deg = width
p.relaxation_s = relax
p.angular_diffusion_s = diffusion
p.radius_mm = radius
p.mass_g = mass
p.drag_kg_s = drag * 1e-3

st.write(f'**Coeficiente de espalhamento:** S = γwa − γia − γiw = **{st.session_state.engine.spreading_coefficient_mNm:+.1f} mN/m**.')
st.caption('O coeficiente S não fixa o gradiente de tensão. O contraste Δγ é um parâmetro efetivo independente.')

@st.fragment(run_every='0.3s')
def live():
    sim = st.session_state.engine
    if st.session_state.running:
        sim.step(dt)
    fx, fy = sim.net_force()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Tempo', f'{sim.t:.2f} s')
    c2.metric('Força longitudinal', f'{fx * 1e6:.2f} μN')
    c3.metric('Velocidade', f'{sim.v * 1e3:.2f} mm/s')
    c4.metric('Deslocamento', f'{sim.x * 1e3:.2f} mm')
    left, right = st.columns(2)
    with left:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(np.rad2deg(sim.theta), sim.gamma_profile() * 1e3)
        ax.set(xlabel='Ângulo na borda (graus)', ylabel='γ efetiva água–ar (mN/m)',
               title='Tensão superficial ao redor do barco', xlim=(0, 360),
               ylim=(max(0, wa - min(dg, wa) - 2), wa + 2))
        ax.grid(alpha=0.3)
        st.pyplot(fig)
        plt.close(fig)
    with right:
        fig, ax = plt.subplots(figsize=(5, 4))
        if sim.history:
            hist = np.asarray(sim.history)
            ax.plot(hist[:, 0], hist[:, 2] * 1e3)
        ax.set(xlabel='Tempo (s)', ylabel='Velocidade (mm/s)', title='Resposta dinâmica do barco')
        ax.grid(alpha=0.3)
        st.pyplot(fig)
        plt.close(fig)
    st.caption('Força capilar integrada na borda circular; movimento longitudinal com arrasto linear.')

live()
