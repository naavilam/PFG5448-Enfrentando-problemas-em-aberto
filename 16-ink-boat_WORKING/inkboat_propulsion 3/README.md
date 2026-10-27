# Ink Boat — simulador de propulsão por tensão superficial

Protótipo Streamlit de **modelo reduzido** para força capilar e movimento longitudinal de um barco circular. Não resolve Navier–Stokes multifásico, a linha de contato trifásica ou o espalhamento da tinta.

## Executar

```bash
pip install -r requirements.txt
streamlit run app.py
```

Para publicar no Streamlit Community Cloud, aponte o *main file path* para `inkboat_propulsion/app.py` (se essa for a pasta no repositório).

## Interface

- **Esquerda:** controles das tensões interfaciais, fonte angular e propriedades do barco.
- **Centro:** força, velocidade, deslocamento e gráficos atualizados.
- **Direita:** todas as equações efetivamente utilizadas pelo modelo, com parâmetros e limitações.

## Equações

A concentração/contraste angular adimensional `c(θ,t)` obedece a

`∂t c = k q(θ) (1-c) - c/τ + Dθ ∂²θ c`,

com `q(θ) = exp[-d₂π(θ,θ₀)²/(2σθ²)]`, e a tensão efetiva é `γ(θ,t)=γwa-Δγ c(θ,t)`.

A força é `F = R ∫ γ(θ,t)(cos θ, sin θ) dθ` e o movimento é `m dv/dt = Fx - b v`, `dx/dt=v`.

O coeficiente de espalhamento `S=γwa-γia-γiw` é exibido, mas não determina `Δγ`. A distribuição de tensão é fenomenológica e necessita calibração experimental.
