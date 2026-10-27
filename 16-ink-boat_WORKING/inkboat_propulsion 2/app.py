"""Interactive reduced Ink Boat propulsion model with governing-equations panel."""
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from model import LiveSimulation

st.set_page_config(page_title='Ink Boat — Propulsão Marangoni', layout='wide')
st.title('Ink Boat | Propulsão por tensão superficial')
st.info('Modelo reduzido: calcula a força capilar resultante na borda e a velocidade do barco. Não resolve o escoamento multifásico nem o espalhamento da tinta.')

with st.sidebar:
    st.header('Interfaces (mN/m)')
    wa = st.slider('Água–ar γwa', 20., 80., 72., 1.)
    ia = st.slider('Tinta–ar γia', 10., 70., 30., 1.)
    iw = st.slider('Tinta–água γiw', 1., 60., 25., 1.)
    st.header('Contraste efetivo na borda')
    dg = st.slider('Redução local Δγ (mN/m)', 0., 40., 8., 0.5)
    release = st.slider('Intensidade da fonte (1/s)', 0., 5., 0.8, 0.1)
    angle = st.slider('Ângulo da fonte (graus)', 0, 359, 180)
    width = st.slider('Largura angular (graus)', 5, 180, 35)
    relax = st.slider('Tempo de relaxação (s)', 0.1, 20., 2., 0.1)
    diffusion = st.slider('Mistura angular (1/s)', 0., 1., 0.15, 0.01)
    st.header('Barco e resistência')
    radius = st.slider('Raio (mm)', 1., 15., 4., 0.5)
    mass = st.slider('Massa (g)', 0.05, 5., 0.5, 0.05)
    drag = st.slider('Arrasto linear b (g/s)', 0.1, 30., 4., 0.1)
    dt = st.slider('Avanço por atualização (s)', 0.01, 0.20, 0.04, 0.01)
    col_play, col_pause, col_reset = st.columns(3)
    play = col_play.button('▶', use_container_width=True)
    pause = col_pause.button('⏸', use_container_width=True)
    reset = col_reset.button('↺', use_container_width=True)

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

st.write(f'**Coeficiente de espalhamento termodinâmico:** '
         f'S = γwa − γia − γiw = **{st.session_state.engine.spreading_coefficient_mNm:+.1f} mN/m**.')
st.caption('As três tensões são parâmetros distintos. S indica uma tendência termodinâmica de espalhamento, '
           'mas não determina a força propulsiva. O contraste Δγ é independente e precisa ser calibrado.')

main, equations = st.columns([3.5, 1.4], gap='large')

with equations:
    st.subheader('Equações do modelo')
    st.caption('Expressões implementadas em model.py; hipóteses de um modelo efetivo, não de Navier–Stokes 3D.')
    st.markdown('**1. Três interfaces**')
    st.latex(r'\gamma_{wa},\quad\gamma_{ia},\quad\gamma_{iw}')
    st.latex(r'S=\gamma_{wa}-\gamma_{ia}-\gamma_{iw}')
    st.markdown('**2. Fonte angular prescrita**')
    st.latex(r'q(\theta)=\exp\!\left[-\frac{d_{2\pi}(\theta,\theta_0)^2}{2\sigma_\theta^2}\right]')
    st.caption('d₂π é a menor distância angular periódica; θ₀ é a direção da fonte.')
    st.markdown('**3. Evolução do contraste superficial**')
    st.latex(r'\frac{\partial c}{\partial t}=k\,q(\theta)(1-c)-\frac{c}{\tau}+D_\theta\frac{\partial^2 c}{\partial\theta^2}')
    st.caption('c é adimensional (0 ≤ c ≤ 1); a mistura é angular e fenomenológica, não uma difusividade molecular.')
    st.markdown('**4. Tensão superficial efetiva**')
    st.latex(r'\gamma(\theta,t)=\gamma_{wa}-\Delta\gamma\,c(\theta,t)')
    st.markdown('**5. Força capilar resultante**')
    st.latex(r'\mathbf F=R\int_0^{2\pi}\gamma(\theta,t)\begin{pmatrix}\cos\theta\\\sin\theta\end{pmatrix}\,d\theta')
    st.latex(r'F_x=R\int_0^{2\pi}\gamma(\theta,t)\cos\theta\,d\theta')
    st.markdown('**6. Movimento longitudinal**')
    st.latex(r'm\,\frac{dv}{dt}=F_x-bv')
    st.latex(r'\frac{dx}{dt}=v')
    with st.expander('Parâmetros e unidades', expanded=False):
        st.markdown(
            f'- $R$ = **{radius:.1f} mm** · $m$ = **{mass:.2f} g**  \n'
            f'- $b$ = **{drag:.1f} g/s** · $\\Delta\\gamma$ = **{dg:.1f} mN/m**  \n'
            f'- $k$ = **{release:.2f} s⁻¹** · $\\tau$ = **{relax:.1f} s**  \n'
            f'- $D_\\theta$ = **{diffusion:.2f} s⁻¹** · $\\theta_0$ = **{angle}°**  \n'
            f'- $\\sigma_\\theta$ = **{width}°**'
        )
    st.warning('A distribuição de tensão é prescrita por uma fonte angular. '
               'As três tensões interfaciais não determinam sozinhas essa distribuição; '
               'o modelo não resolve a linha de contato trifásica.')

with main:
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
            theta = sim.theta
            gamma = sim.gamma_profile() * 1e3
            ax.plot(np.rad2deg(theta), gamma)
            ax.set(xlabel='Ângulo na borda (graus)',
                   ylabel='γ efetiva água–ar (mN/m)',
                   title='Tensão superficial ao redor do barco',
                   xlim=(0, 360), ylim=(max(0, wa - dg - 2), wa + 2))
            ax.grid(alpha=0.3)
            st.pyplot(fig)
            plt.close(fig)
        with right:
            fig, ax = plt.subplots(figsize=(5, 4))
            if sim.history:
                hist = np.asarray(sim.history)
                ax.plot(hist[:, 0], hist[:, 2] * 1e3)
            ax.set(xlabel='Tempo (s)', ylabel='Velocidade (mm/s)',
                   title='Resposta dinâmica do barco')
            ax.grid(alpha=0.3)
            st.pyplot(fig)
            plt.close(fig)
        st.caption('Força capilar de uma linha de contato circular; movimento apenas em x; arrasto linear.')

    live()
