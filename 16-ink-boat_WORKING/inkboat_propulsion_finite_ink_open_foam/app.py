"""Ink boat — Marangoni 3D (OpenFOAM), interactive viewer of a precomputed run.

The solver works in the boat frame; here everything is shown in the LAB frame:
x_lab = x + X_boat(t),  u_lab = u + U_boat(t) e_x.
"""
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st

DATA = Path(__file__).with_name("data") / "marangoni_fields.npz"
INK = "Blues"
WOOD = "#8d5a2b"

st.set_page_config(page_title="Ink boat · Marangoni 3D", layout="wide")


@st.cache_data
def load():
    d = np.load(DATA)
    out = {k: d[k] for k in d.files}
    for k, v in out.items():
        if v.dtype == np.float16:
            out[k] = v.astype(np.float32)
    return out


D = load()
x, y, z, t = D["x"], D["y"], D["z"], D["t"]
bx0, bx1, by0, by1, bz0, bz1 = D["boat"]
Xf, Uf = D["frame_X"], D["frame_U"]          # boat position / velocity at each frame
cm = 100.0
ink_end = float(D["ink_duration"])
x_range = [(x[0]) * cm, (x[-1] + Xf.max()) * cm]


# ----------------------------------------------------------------------------- helpers
def boat_mesh(shift, z_scale=1.0):
    X = (np.array([bx0, bx1, bx1, bx0, bx0, bx1, bx1, bx0]) + shift) * cm
    Y = np.array([by0, by0, by1, by1, by0, by0, by1, by1]) * cm
    top = bz1 + 0.0004
    Z = np.array([bz0, bz0, bz0, bz0, top, top, top, top]) * cm * z_scale
    i = [0, 0, 4, 4, 0, 0, 1, 1, 2, 2, 3, 3]
    j = [1, 2, 5, 6, 1, 5, 2, 6, 3, 7, 0, 4]
    k = [2, 3, 6, 7, 5, 4, 6, 5, 7, 6, 4, 7]
    return go.Mesh3d(x=X, y=Y, z=Z, i=i, j=j, k=k, color=WOOD, opacity=1.0, flatshading=True,
                     name="barco", hoverinfo="name", showscale=False)


def start_marker(z_scale=1.0):
    """Fixed post at the boat's initial bow position: makes the motion visible."""
    return go.Scatter3d(x=[bx1 * cm, bx1 * cm], y=[y[0] * cm, y[0] * cm], z=[0, 1.2 * z_scale],
                        mode="lines+text", text=["", "início"], textposition="top center",
                        line=dict(color="#555", width=6), hoverinfo="skip", showlegend=False)


def arrows_2d(X, Y, U, V, scale, min_frac=0.03):
    xs, ys = [], []
    sp = np.hypot(U, V)
    smax = np.nanmax(sp) if np.isfinite(sp).any() else 0
    for x0, y0, u, v, s in zip(X.ravel(), Y.ravel(), U.ravel(), V.ravel(), sp.ravel()):
        if not np.isfinite(s) or s <= min_frac * smax:
            continue
        x1, y1 = x0 + scale * u, y0 + scale * v
        ang = np.arctan2(v, u)
        h = 0.3 * scale * s
        xs += [x0, x1, None, x1 - h * np.cos(ang - 0.4), x1, x1 - h * np.cos(ang + 0.4), None]
        ys += [y0, y1, None, y1 - h * np.sin(ang - 0.4), y1, y1 - h * np.sin(ang + 0.4), None]
    return xs, ys


def frame_label(i):
    return f"{t[i]:.1f} s"


def anim_controls(n, slider_y, button_y):
    return dict(
        updatemenus=[dict(type="buttons", x=0.0, y=button_y, xanchor="left", buttons=[
            dict(label="▶", method="animate",
                 args=[None, dict(frame=dict(duration=150, redraw=True), fromcurrent=True, transition=dict(duration=0))]),
            dict(label="❚❚", method="animate", args=[[None], dict(mode="immediate", frame=dict(duration=0))]),
        ])],
        sliders=[dict(x=0.1, len=0.9, y=slider_y, currentvalue=dict(prefix="t = "), steps=[
            dict(label=frame_label(i), method="animate",
                 args=[[frame_label(i)], dict(mode="immediate", frame=dict(duration=0, redraw=True))])
            for i in range(n)])],
    )


# ----------------------------------------------------------------------------- header
st.title("Ink boat — propulsão Marangoni 3D (OpenFOAM)")
st.caption(
    f"Uma gota finita de tinta é liberada na popa por {ink_end:.0f} s. A tinta reduz a tensão superficial; "
    "a diferença γ_frente − γ_trás puxa o barco para a frente e o gradiente de γ arrasta a água para longe da popa. "
    "A força resultante (linha de contato + pressão + viscosa) acelera o barco, m dU/dt = Fₓ. "
    "Simulação pré-calculada, mostrada no referencial do laboratório: gire, aproxime e anime."
)

c1, c2, c3, c4 = st.columns(4)
c1.metric("Velocidade máxima do barco", f"{100*D['U_boat'].max():.2f} cm/s")
c2.metric("Distância percorrida", f"{100*D['X_boat'][-1]:.1f} cm")
c3.metric("Força de linha de contato máx.", f"{1e6*D['F_cl_x'].max():.1f} µN")
c4.metric("Massa do barco", f"{1e3*float(D['boat_mass']):.2f} g")

tab_surf, tab_motion, tab_vol, tab_slice, tab_model = st.tabs(
    ["Superfície 3D", "Movimento e forças", "Volume 3D", "Corte vertical", "Modelo"])

# ----------------------------------------------------------------------------- surface (animated, lab frame)
with tab_surf:
    st.markdown("**Tinta na superfície e velocidade da água (referencial do laboratório).** "
                "▶ anima; o poste cinza marca a posição inicial da proa. Gire a cena com o mouse.")
    Xs, Ys = np.meshgrid(x * cm, y * cm, indexing="ij")
    sub = (slice(2, None, 5), slice(2, None, 5))
    u_lab = D["surf_ux"] + Uf[:, None, None]
    vmax = float(np.nanmax(np.hypot(u_lab, D["surf_uy"])))
    a_scale = 1.2 / max(vmax, 1e-9)
    zs = 3.0

    def surf_traces(i):
        shift = Xf[i] * cm
        ax, ay = arrows_2d(Xs[sub] + shift, Ys[sub], u_lab[i][sub], D["surf_uy"][i][sub], a_scale)
        return [
            go.Surface(x=Xs + shift, y=Ys, z=np.zeros_like(Xs), surfacecolor=np.nan_to_num(D["surf_c"][i], nan=0.0),
                       cmin=0, cmax=1, colorscale=INK, colorbar=dict(title="c", len=0.6),
                       hovertemplate="x=%{x:.1f} cm<br>y=%{y:.1f} cm<br>c=%{surfacecolor:.2f}<extra></extra>"),
            go.Scatter3d(x=ax, y=ay, z=[0.03 if v is not None else None for v in ax], mode="lines",
                         line=dict(color="#d1495b", width=3), hoverinfo="skip"),
            boat_mesh(Xf[i], z_scale=zs),
        ]

    fig = go.Figure(data=surf_traces(0) + [start_marker(zs)],
                    frames=[go.Frame(data=surf_traces(i), traces=[0, 1, 2], name=frame_label(i))
                            for i in range(len(t))])
    fig.update_layout(
        height=620, margin=dict(l=0, r=0, t=10, b=0), showlegend=False,
        scene=dict(xaxis=dict(title="x (cm)", range=x_range), yaxis_title="y (cm)",
                   zaxis=dict(title="z (cm, ×3)", range=[z[0] * cm * zs, 2.0 * zs]),
                   aspectmode="manual",
                   aspectratio=dict(x=(x_range[1] - x_range[0]) / 16, y=0.5, z=0.25),
                   camera=dict(eye=dict(x=-0.4, y=-1.6, z=0.9))),
        **anim_controls(len(t), slider_y=0.05, button_y=0.1),
    )
    st.plotly_chart(fig, width="stretch")

# ----------------------------------------------------------------------------- motion
with tab_motion:
    ft = D["force_t"]
    col_a, col_b = st.columns(2)
    with col_a:
        fig = go.Figure([
            go.Scatter(x=ft, y=100 * D["U_boat"], name="U (cm/s)", line=dict(color="#222", width=3)),
            go.Scatter(x=ft, y=100 * D["X_boat"], name="X (cm)", line=dict(color="#2e86ab", width=2.5, dash="dash"),
                       yaxis="y2"),
        ])
        fig.add_vline(x=ink_end, line=dict(color="#999", dash="dot"), annotation_text="tinta acaba")
        fig.update_layout(height=400, title="Movimento do barco", xaxis_title="t (s)",
                          yaxis=dict(title="velocidade U (cm/s)"),
                          yaxis2=dict(title="posição X (cm)", overlaying="y", side="right"),
                          legend=dict(orientation="h", y=-0.2), margin=dict(l=10, r=10, t=40, b=10))
        st.plotly_chart(fig, width="stretch")
    with col_b:
        fig = go.Figure([
            go.Scatter(x=ft, y=1e6 * D["F_cl_x"], name="linha de contato ∮γ t dℓ", line=dict(color="#2e86ab", width=2.5)),
            go.Scatter(x=ft, y=1e6 * D["F_pressure_x"], name="pressão", line=dict(color="#edae49", width=2)),
            go.Scatter(x=ft, y=1e6 * D["F_viscous_x"], name="viscosa", line=dict(color="#d1495b", width=2)),
            go.Scatter(x=ft, y=1e6 * D["F_net_x"], name="líquida = m dU/dt", line=dict(color="#222", width=3, dash="dash")),
        ])
        fig.add_vline(x=ink_end, line=dict(color="#999", dash="dot"))
        fig.add_hline(y=0, line=dict(color="#ccc"))
        fig.update_layout(height=400, title="Forças sobre o barco (+ = para a frente)", xaxis_title="t (s)",
                          yaxis_title="Fₓ (µN)", legend=dict(orientation="h", y=-0.2),
                          margin=dict(l=10, r=10, t=40, b=10))
        st.plotly_chart(fig, width="stretch")
    st.markdown(
        "Balanço da Eq. 3.5 da tese: **m dU/dt = ∮ γ t dℓ + ∫ n·σ dS**. "
        "Com a tinta, a linha de contato empurra e o casco é freado por pressão e cisalhamento; "
        "quando a tinta acaba, o contraste de γ some e o arrasto desacelera o barco.")

# ----------------------------------------------------------------------------- volume
with tab_vol:
    st.markdown("**Pluma de tinta em 3D** (isosuperfícies de c) e velocidade da água no referencial do "
                "laboratório. Escala vertical ampliada 3×.")
    col_a, col_b = st.columns([3, 1])
    with col_b:
        ti = st.select_slider("Tempo", options=list(range(len(t))), value=min(12, len(t) - 1),
                              format_func=frame_label, key="vol_t")
        show_cones = st.toggle("Mostrar velocidade", value=True)
        levels = st.select_slider("Faixa de c", options=[0.02, 0.05, 0.1, 0.2, 0.4], value=(0.05, 0.4))
    zs3 = 3.0
    shift = Xf[ti]
    X3, Y3, Z3 = np.meshgrid((D["x2"] + shift) * cm, D["y2"] * cm, D["z2"] * cm * zs3, indexing="ij")
    vc = np.nan_to_num(D["vol_c"][ti], nan=0.0)
    traces = [go.Isosurface(x=X3.ravel(), y=Y3.ravel(), z=Z3.ravel(), value=vc.ravel(),
                            isomin=levels[0], isomax=levels[1], surface_count=3, opacity=0.35,
                            colorscale=INK, cmin=0, cmax=0.6, caps=dict(x_show=False, y_show=False, z_show=False),
                            colorbar=dict(title="c", len=0.6)),
              boat_mesh(shift, z_scale=zs3), start_marker(zs3)]
    if show_cones:
        ux = np.nan_to_num(D["vol_ux"][ti], nan=0.0) + Uf[ti]
        uy = np.nan_to_num(D["vol_uy"][ti], nan=0.0)
        uz = np.nan_to_num(D["vol_uz"][ti], nan=0.0)
        boat_cells = np.isnan(D["vol_ux"][ti])
        ux[boat_cells] = 0.0
        sl = (slice(1, None, 3), slice(1, None, 3), slice(None, None, 2))
        sp = np.sqrt(ux**2 + uy**2 + uz**2)[sl]
        keep = sp.ravel() > 0.08 * max(sp.max(), 1e-12)
        traces.append(go.Cone(x=X3[sl].ravel()[keep], y=Y3[sl].ravel()[keep], z=Z3[sl].ravel()[keep],
                              u=ux[sl].ravel()[keep], v=uy[sl].ravel()[keep], w=zs3 * uz[sl].ravel()[keep],
                              sizemode="scaled", sizeref=1.2, anchor="tail", colorscale="Reds",
                              showscale=False, opacity=0.8, hoverinfo="skip"))
    fig = go.Figure(traces)
    fig.update_layout(height=620, margin=dict(l=0, r=0, t=10, b=0),
                      scene=dict(xaxis=dict(title="x (cm)", range=x_range), yaxis_title="y (cm)",
                                 zaxis_title="z (cm, ×3)", aspectmode="manual",
                                 aspectratio=dict(x=(x_range[1] - x_range[0]) / 16, y=0.5, z=0.4),
                                 camera=dict(eye=dict(x=-0.5, y=-1.7, z=0.8))))
    with col_a:
        st.plotly_chart(fig, width="stretch")

# ----------------------------------------------------------------------------- vertical slice (animated, lab frame)
with tab_slice:
    st.markdown("**Plano central y = 0, referencial do laboratório.** A tensão Marangoni puxa a superfície "
                "para trás da popa, e a continuidade fecha o circuito com um retorno mais fundo. "
                "Eixo vertical ampliado.")
    Xm, Zm = np.meshgrid(x * cm, z * cm, indexing="ij")
    subm = (slice(1, None, 4), slice(None, None, 2))
    ux_m = D["slice_ux"] + Uf[:, None, None]
    vmax_m = float(np.nanmax(np.hypot(ux_m, D["slice_uz"])))
    sc_m = 1.0 / max(vmax_m, 1e-9)

    def slice_traces(i):
        shift = Xf[i] * cm
        ax, az = arrows_2d(Xm[subm] + shift, Zm[subm], ux_m[i][subm], D["slice_uz"][i][subm], sc_m)
        bx = np.array([bx0, bx1, bx1, bx0, bx0]) * cm + shift
        bz = np.array([bz0, bz0, 0, 0, bz0]) * cm
        return [go.Heatmap(x=x * cm + shift, y=z * cm, z=D["slice_c"][i].T, zmin=0, zmax=1, colorscale=INK,
                           colorbar=dict(title="c"),
                           hovertemplate="x=%{x:.1f} cm<br>z=%{y:.2f} cm<br>c=%{z:.2f}<extra></extra>"),
                go.Scatter(x=ax, y=az, mode="lines", line=dict(color="#d1495b", width=1.4), hoverinfo="skip"),
                go.Scatter(x=bx, y=bz, fill="toself", fillcolor=WOOD, line=dict(color=WOOD), hoverinfo="skip")]

    fig = go.Figure(data=slice_traces(0),
                    frames=[go.Frame(data=slice_traces(i), name=frame_label(i)) for i in range(len(t))])
    fig.add_vline(x=bx1 * cm, line=dict(color="#555", dash="dot"), annotation_text="início")
    fig.update_layout(height=430, margin=dict(l=10, r=10, t=10, b=10), showlegend=False,
                      xaxis=dict(title="x (cm)", range=x_range), yaxis=dict(title="z (cm)", range=[z[0] * cm, 0.05]),
                      **anim_controls(len(t), slider_y=-0.12, button_y=-0.18))
    st.plotly_chart(fig, width="stretch")

# ----------------------------------------------------------------------------- model
with tab_model:
    st.markdown(rf"""
**Solver:** OpenFOAM v2606, `pimpleFoam` (Navier–Stokes incompressível, laminar), transiente, 8 processos.

**Barco livre em surge (translação em x):** $m\,dU/dt = F_x$, com
$F_x = \oint\gamma\,\mathbf t\,d\ell + \int \mathbf n\cdot\boldsymbol\sigma\,dS$ (Eq. 3.5 da tese).
$m = {1e3*float(D['boat_mass']):.2f}$ g, igual à água deslocada por um casco de 3 × 1,5 cm e 3 mm de calado.
A integração é explícita, a cada passo de tempo, num *function object* `coded`.

**Referencial do barco:** a malha (janela de 16 × 8 × 2 cm) acompanha o barco. A água entra pela proa a $-U(t)$,
o fundo da cuba se move a $-U(t)$ e a aceleração do referencial entra como força fictícia $-a$ na água.
Este app converte tudo para o laboratório: $x_\text{{lab}} = x + X(t)$, $u_\text{{lab}} = u + U(t)$.

**Tinta:** escalar passivo $c\in[0,1]$ (advecção–difusão, $D = 2\times10^{{-6}}$ m²/s efetiva).
A fonte presa à popa relaxa $c\to1$ durante {ink_end:.0f} s e então desliga (**tinta finita**).

**Marangoni:** $\gamma = \gamma_0 - \Gamma c$ com $\Gamma = {1e3*float(D['gamma_coeff']):.1f}$ mN/m. A superfície livre é plana
(Ca ≈ 10⁻⁴), com deslizamento. A tensão tangencial $\nabla_s\gamma$ entra como força de corpo na camada de células
da superfície. A força de linha de contato soma $-\Gamma c\,\mathbf t\,d\ell$ nas células vizinhas ao casco;
$\gamma_0$ uniforme integra a zero.

**Limites:** apenas surge (sem guinada nem arfagem); superfície não deforma; a tinta é um escalar de volume, não
um surfactante só de superfície com adsorção; γ(c) linear; o cisalhamento no casco omite a parte normal do
tensor viscoso (diferença de ~2% em relação ao `forces` do OpenFOAM). Os valores são ilustrativos.
""")
