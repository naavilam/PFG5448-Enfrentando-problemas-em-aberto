"""Run with: python -m streamlit run app.py"""
import csv
import io
import time
from dataclasses import astuple, replace
from playback import Playback

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from model import (
    Params, WATER_DENSITY_KG_M3, WATER_PRESSURE_PA, WATER_SOURCE,
    WATER_TEMPERATURE_C, WATER_VISCOSITY_PA_S, contact_line, gamma_on_contour, simulate,
)

st.set_page_config(page_title="Ink Boat | Balanço de forças", layout="wide")
st.markdown(
    "<style>[data-testid=stSidebar]{min-width:610px!important;width:610px!important;}"
    "@media(max-width:900px){[data-testid=stSidebar]{min-width:350px!important;width:350px!important;}}"
    "</style>", unsafe_allow_html=True,
)
st.title("Ink Boat — Balanço de forças")
st.markdown("**Equação 3.5 · contraste de tensão superficial e resposta viscosa de um disco fino**")
st.caption("Equações, constantes e controles à esquerda; velocidade, aceleração e forças à direita.")

with st.sidebar:
    st.header("Simulação")
    c1, c2, c3 = st.columns(3)
    start = c1.button("▶ Iniciar / Retomar", width="stretch")
    pause = c2.button("⏸ Pausar", width="stretch")
    stop = c3.button("⏹ Parar / Reiniciar", width="stretch")
    playback_rate = st.select_slider("Velocidade de reprodução", options=[0.1, 0.25, 0.5, 1.0, 2.0, 5.0], value=1.0,
                                     format_func=lambda value: f"{value:g}×")
    st.caption("Pausar mantém o instante atual. Parar retorna a t=0. Alterar parâmetros físicos reinicia o ensaio.")
    st.header("Equações e parâmetros")
    with st.container(border=True):
        st.markdown("**1 · Balanço de quantidade de movimento — Eq. 3.5**")
        st.latex(r"\oint_{\ell_p}\gamma\,\mathbf t\,d\ell+\int_{S_p}\mathbf n\cdot\boldsymbol\sigma\,dS=m\frac{d\mathbf U}{dt}")
        st.latex(r"\mathbf U(0)=\mathbf0")
        st.caption("Translação horizontal em x. A força de peso e o empuxo não entram nessa componente; a imersão é prescrita.")
        mass_g = st.slider("Massa do barco m (g)", 0.02, 5.0, 0.20, 0.02)
        st.latex(r"m=" + f"{mass_g*1e-3:.3g}" + r"\;\mathrm{kg}")

    with st.container(border=True):
        st.markdown("**2 · Contorno do barco e tensão superficial**")
        st.latex(r"P=\oint_{\ell_p}d\ell=2\pi R,\qquad A_w=\pi R^2")
        perimeter_mm = st.slider("Perímetro total P (mm)", 2.0, 100.0, 20.0, 1.0)
        gamma_mnm = st.slider("Tensão superficial média γ₀ (mN/m)", 0.0, 100.0, 72.0, 0.5)
        contrast_limit = max(2*gamma_mnm, 0.001)
        contrast_mnm = st.number_input("Contraste inicial Δγ₀ (mN/m)",
                                       min_value=-contrast_limit, max_value=contrast_limit,
                                       value=min(0.01, 2*gamma_mnm), step=0.001,
                                       format="%.4f", disabled=gamma_mnm == 0)
        if gamma_mnm == 0:
            contrast_mnm = 0.0
        st.latex(r"\gamma(\theta,t)=\gamma_0+\frac{\Delta\gamma(t)}{2}\cos\theta")
        st.latex(r"\gamma_{\rm frente}=\gamma_0+\frac{\Delta\gamma(t)}{2},\quad\gamma_{\rm tr\acute as}=\gamma_0-\frac{\Delta\gamma(t)}{2}")
        st.caption(f"Frente inicial (+x): {gamma_mnm+contrast_mnm/2:.4f} mN/m · Traseira inicial (−x): {gamma_mnm-contrast_mnm/2:.4f} mN/m")
        st.latex(r"F_{\gamma,x}=R\int_0^{2\pi}\gamma(\theta,t)\cos\theta\,d\theta=\frac{P\Delta\gamma(t)}{4}")
        st.caption("ℓₚ é a linha de contato; P é seu comprimento. Na tese, t é tangente à interface água–ar e normal ao contorno. Para o disco, aponta radialmente para fora.")
        st.caption("Δγ₀ > 0 representa menor tensão atrás e força para +x. A média γ₀ fica fixa; o contraste decai com a fração de tinta restante, até desaparecer.")

    with st.container(border=True):
        st.markdown("**3 · Gota de tinta e contraste no tempo**")
        volume_ul = st.number_input("Volume inicial da gota V₀ (µL)", min_value=0.0, max_value=1000.0,
                                    value=0.02, step=0.005, format="%.4f")
        rate_ul_s = st.number_input("Taxa de ejeção Q (µL/s)", min_value=0.0, max_value=100.0,
                                    value=0.005, step=0.001, format="%.4f")
        st.latex(r"V(t)=\max(V_0-Qt,0),\quad t_{\rm fim}=V_0/Q")
        st.latex(r"\Delta\gamma(t)=\Delta\gamma_0\frac{V(t)}{V_0}")
        if volume_ul == 0:
            st.caption("Sem tinta: contraste e força propulsiva nulos desde t=0.")
        elif rate_ul_s == 0:
            st.caption("Q=0: sem consumo; o contraste inicial permanece constante.")
        else:
            st.caption(f"Tinta esgotada em {volume_ul/rate_ul_s:.3g} s. Depois disso, a força capilar é zero e o barco desacelera pelo arrasto.")
        st.caption("Hipótese inicial: contraste proporcional ao volume restante. Não é uma lei de transporte/difusão da tinta; γ₀ não vai a zero.")

    with st.container(border=True):
        st.markdown("**4 · Tensor de tensões — Eq. 3.3**")
        st.latex(r"\boldsymbol\sigma=-p\mathbf I+\mu\left[\nabla\mathbf u+(\nabla\mathbf u)^T\right]")
        st.latex(r"\mu=" + f"{WATER_VISCOSITY_PA_S:.7f}" + r"\;\mathrm{Pa\,s}")
        st.markdown(f"**Constante:** μ = {WATER_VISCOSITY_PA_S*1e3:.4f} mPa·s, água pura a {WATER_TEMPERATURE_C:.0f} °C e ≈ {WATER_PRESSURE_PA/1e6:.1f} MPa.")
        st.caption("ρ ≈ 998,2 kg/m³ é usada apenas no diagnóstico de Reynolds. As propriedades da água permanecem fixas nesta versão.")
        st.markdown(f"[Fonte: IAPWS SR6-08(2011), seção 3, Eq. (7)]({WATER_SOURCE})")
        st.markdown("**Fechamento adotado: cisalhamento sob a face inferior**")
        st.latex(r"u_x(z,t)\simeq U(t)\left(1+\frac zh\right),\quad-h\le z\le0")
        st.latex(r"F_{h,x}\simeq-\mu\frac{A_w}{h}U=-bU,\qquad b=\frac{\mu A_w}{h}")
        shear_depth_mm = st.slider("Espessura efetiva de cisalhamento h (mm)", 0.05, 5.0, 0.50, 0.05)
        st.caption("Água em repouso em z=−h; não deslizamento na face z=0. Pressão sem gradiente horizontal, sem contribuição horizontal líquida. Área lateral do disco desprezada.")
        st.caption("h é uma escala efetiva do fechamento, não uma profundidade medida. Não estamos calculando os campos p(x,y,z,t) e u(x,y,z,t) da tese.")

    with st.container(border=True):
        st.markdown("**5 · Velocidade e aceleração no tempo**")
        st.latex(r"m\dot U=F_\gamma-bU,\qquad U(0)=0")
        st.latex(r"\tau=\frac mb,\qquad a(t)=\frac{F_\gamma(t)-bU(t)}m")
        st.latex(r"t>t_{\rm fim}:\quad U(t)=U(t_{\rm fim})e^{-(t-t_{\rm fim})/\tau}")
        duration_s = st.slider("Janela de tempo (s)", 0.1, 60.0, 5.0, 0.1)
        points = st.select_slider("Pontos nas curvas", options=[101, 251, 501, 1001, 2001], value=501)
        st.caption("Solução analítica por trechos: força decrescente durante o consumo e relaxação viscosa após o esgotamento.")

p = Params(mass_kg=mass_g*1e-3, perimeter_m=perimeter_mm*1e-3,
           gamma_n_m=gamma_mnm*1e-3, gamma_contrast_n_m=contrast_mnm*1e-3, shear_depth_m=shear_depth_mm*1e-3,
           duration_s=duration_s, points=points,
           ink_volume_m3=volume_ul*1e-9, ejection_rate_m3_s=rate_ul_s*1e-9)
signature = astuple(p)
now = time.monotonic()
if st.session_state.get("inkboat_parameters") != signature:
    st.session_state.inkboat_clock = Playback(last_tick=now)
    st.session_state.inkboat_parameters = signature
clock = st.session_state.inkboat_clock
clock.advance(now, p.duration_s, playback_rate)
if stop:
    clock.stop(now)
elif pause:
    clock.pause()
elif start:
    clock.start(now, p.duration_s)


if p.gamma_contrast_n_m == 0 or p.ink_volume_m3 == 0:
    st.info("Contraste zero: γ uniforme ⇒ força capilar nula. Partindo do repouso, U(t)=0 e a(t)=0.")
else:
    st.info("Fγ(t) = P Δγ(t) / 4. O contraste diminui com o consumo da tinta; "
            "quando a tinta acaba, a propulsão cessa e a velocidade decai por arrasto.")

@st.fragment(run_every=0.25)
def live():
    clock = st.session_state.inkboat_clock
    clock.advance(time.monotonic(), p.duration_s, playback_rate)
    elapsed = clock.elapsed_s
    if elapsed > 0:
        count = max(2, int((p.points-1)*elapsed/p.duration_s)+1)
        r = simulate(replace(p, duration_s=elapsed, points=count))
    else:
        r = {key: values[:1] for key, values in simulate(p).items()}
    status = "Em execução" if clock.running else ("Concluída" if elapsed >= p.duration_s else ("Pausada" if elapsed > 0 else "Pronta"))
    st.markdown(f"**{status}** · Tempo simulado: **{elapsed:.2f} / {p.duration_s:.2f} s**")
    st.progress(min(1.0, elapsed/p.duration_s))
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Velocidade atual U", f"{r['velocity_m_s'][-1]*1e3:.4f} mm/s")
    m2.metric("Aceleração atual", f"{r['acceleration_m_s2'][-1]*1e3:.4f} mm/s²")
    m3.metric("Força capilar resultante", f"{r['capillary_force_n'][-1]*1e6:.4f} µN")
    m4.metric("Tempo de relaxação τ", f"{p.relaxation_time_s:.3g} s")
    c1, c2, c3 = st.columns(3)
    c1.metric("Contraste atual Δγ", f"{r['contrast_n_m'][-1]*1e3:.5f} mN/m")
    c2.metric("Tinta restante", f"{r['ink_remaining_m3'][-1]*1e9:.4f} µL")
    c3.metric("Tinta ejetada", f"{r['ink_ejected_m3'][-1]*1e9:.4f} µL")
    if np.isfinite(p.empty_time_s):
        st.caption(f"Esgotamento previsto: {p.empty_time_s:.3g} s · "
                   + ("Reservatório vazio; movimento residual por inércia." if elapsed >= p.empty_time_s else "Reservatório em consumo."))
    st.caption(f"Disco: R={p.radius_m*1e3:.3f} mm · P={perimeter_mm:.2f} mm · "
               f"área inferior molhada={p.wetted_area_m2*1e6:.3f} mm² · b={p.drag_kg_s:.4g} kg/s")


    def curve(key, label, color):
        fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
        ax.plot(r["time_s"], r[key]*1e3, color=color, lw=2.3, marker="o" if len(r["time_s"]) == 1 else None)
        ax.axhline(0, color="#808080", lw=0.6)
        ax.set(xlabel="Tempo (s)", ylabel=label, xlim=(0, p.duration_s))
        if 0 < p.empty_time_s <= p.duration_s:
            ax.axvline(p.empty_time_s, color="gray", linestyle=":", lw=1, label="Fim da tinta")
        ax.grid(alpha=0.22)
        st.pyplot(fig, width="stretch")
        plt.close(fig)

    left, right = st.columns(2, gap="large")
    with left:
        st.subheader("Velocidade do barco")
        curve("velocity_m_s", "Velocidade Uₓ (mm/s)", "#286080")
    with right:
        st.subheader("Aceleração do barco")
        curve("acceleration_m_s2", "Aceleração aₓ (mm/s²)", "#c06c2b")

    left, right = st.columns(2, gap="large")
    with left:
        st.subheader("Tinta e contraste no tempo")
        fig, ax = plt.subplots(figsize=(6.4, 3.2), constrained_layout=True)
        fraction = r["ink_remaining_m3"]/p.ink_volume_m3 if p.ink_volume_m3 > 0 else np.zeros_like(r["time_s"])
        ax.plot(r["time_s"], fraction*100, color="#4b8b53", lw=2)
        ax.set(xlabel="Tempo (s)", ylabel="Tinta restante / contraste inicial (%)", xlim=(0, p.duration_s), ylim=(-2, 105))
        if 0 < p.empty_time_s <= p.duration_s:
            ax.axvline(p.empty_time_s, color="gray", linestyle=":", lw=1)
        ax.grid(alpha=0.22)
        st.pyplot(fig, width="stretch")
        plt.close(fig)
        st.subheader("Forças horizontais")
        fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
        for key, label, style in [
            ("capillary_force_n", "Tensão superficial Fγ", "-"),
            ("hydrodynamic_force_n", "Viscosa Fh", "--"),
            ("total_force_n", "Resultante = m a", "-."),
        ]:
            ax.plot(r["time_s"], r[key]*1e6, style, label=label, lw=2)
        ax.set(xlabel="Tempo (s)", ylabel="Força em x (µN)", xlim=(0, p.duration_s))
        if 0 < p.empty_time_s <= p.duration_s:
            ax.axvline(p.empty_time_s, color="gray", linestyle=":", lw=1, label="Fim da tinta")
        ax.grid(alpha=0.22)
        ax.legend(fontsize=9)
        st.pyplot(fig, width="stretch")
        plt.close(fig)
    with right:
        st.subheader("Contorno: contraste de tensão superficial")
        theta, xy, normals = contact_line(p)
        fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
        closed = np.vstack((xy, xy[0]))*1e3
        ax.fill(closed[:, 0], closed[:, 1], color="#d9e6ec", alpha=0.8)
        ax.plot(closed[:, 0], closed[:, 1], color="#286080", lw=2)
        ids = np.arange(0, len(theta), 16)
        if p.gamma_n_m > 0:
            gamma = gamma_on_contour(theta, elapsed, p)
            arrows = normals[ids]*(gamma[ids]/gamma.max())[:, None]*p.radius_m*0.3*1e3
            ax.quiver(xy[ids, 0]*1e3, xy[ids, 1]*1e3, arrows[:, 0], arrows[:, 1],
                      angles="xy", scale_units="xy", scale=1, color="#c06c2b", width=0.007)
        extent = p.radius_m*1e3*1.6
        ax.set(xlabel="x (mm)", ylabel="y (mm)", aspect="equal",
               xlim=(-extent, extent), ylim=(-extent, extent))
        ax.grid(alpha=0.15)
        st.pyplot(fig, width="stretch")
        plt.close(fig)
        st.caption("Setas indicam a direção local da tração capilar; a escala é ilustrativa. "
                   "A intensidade relativa das setas acompanha γ. A componente uniforme se cancela; o contraste gera a resultante em x.")


    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(r.keys())
    writer.writerows(zip(*r.values()))
    st.download_button("Baixar curvas até o instante atual (CSV, SI)", data=buffer.getvalue(),
                       file_name="inkboat_force_balance.csv", mime="text/csv")

live()

with st.expander("Hipóteses, dados que faltam e próximos passos"):
    st.markdown(r"""
- **Geometria:** disco fino, linha de contato circular plana, área molhada inferior $A_w = \pi R^2$; sem rotação ou variação de imersão.
- **Fechamento viscoso:** perfil linear em uma camada h, pressão sem gradiente horizontal e arrasto $b = \mu A_w/h$. É uma aproximação escolhida, não a integração de um escoamento resolvido.
- **Condição inicial:** barco e água em repouso. O efeito da tinta é representado apenas por Δγ prescrito; há consumo de um reservatório finito, mas sem transporte de concentração.
- **Para o barco real:** faltam formato, massa e imersão medidos, condições nas paredes/fundo, temperatura e uma lei espacial/temporal de γ vinculada à concentração de tinta. Para avaliar a Eq. 3.3 completa, precisamos resolver p e u com essas condições.
- **Próxima etapa:** substituir γ₀ por γ(θ,t). O contraste atual é proporcional ao volume restante e decresce linearmente até o esgotamento. Um gradiente real também pode gerar escoamento de Marangoni na água; esse efeito não é resolvido pelo fechamento de arrasto e exigirá uma revisão na próxima etapa.
""")
    re_max = float(simulate(p)["reynolds"].max())
    diffusion_time = WATER_DENSITY_KG_M3*p.shear_depth_m**2/WATER_VISCOSITY_PA_S
    st.write(f"Diagnósticos da janela completa: Re máximo (diâmetro) = {re_max:.3g}; "
             f"h/R = {p.shear_depth_m/p.radius_m:.3g}; "
             f"tempo de difusão viscosa h²/ν ≈ {diffusion_time:.3g} s.")
    st.caption("Camada fina, efeitos de borda pequenos e difusão viscosa rápida em relação à resposta "
               "do barco favorecem esse fechamento. Não foram incluídos inércia do fluido, massa adicionada ou memória hidrodinâmica.")
    st.markdown("Tese: Saeed Jafari Kang (2021), *Marangoni Propulsion of Active Particles*, "
                "§3.2, Eqs. (3.3) e (3.5), pp. 26–27. PDF preservado nesta pasta.")

