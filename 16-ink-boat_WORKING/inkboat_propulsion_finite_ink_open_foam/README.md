# Ink boat — escoamento Marangoni 3D (OpenFOAM + Streamlit)

Simulação 3D, feita uma única vez offline em OpenFOAM, da propulsão Marangoni de um barco
**livre** (translação em x) por uma gota finita de tinta na popa. A força ∮γ t dℓ + ∫n·σ dS
acelera o barco: m dU/dt = Fₓ. O app Streamlit mostra os
campos pré-calculados de forma interativa: cena 3D girável, animação no tempo,
isosuperfícies, corte vertical e forças. O app não roda CFD e não tem parâmetros físicos.

```
case/                 caso OpenFOAM (v2606, pimpleFoam laminar)
export_fields.py      caso reconstruído  ->  data/marangoni_fields.npz (só numpy)
data/                 campos exportados usados pelo app
app.py                visualizador Streamlit (plotly)
```

## Ver o app

```sh
python3 -m pip install -r requirements.txt
python3 -m streamlit run app.py
```

Para publicar no Streamlit Community Cloud, suba esta pasta (com `data/`) num repositório
GitHub e aponte o app para `app.py`. O NPZ tem poucos MB e não precisa de OpenFOAM no servidor.

## Refazer a simulação

OpenFOAM v2606 (macOS: `brew install --cask gerlero/openfoam/openfoam`). Rode numa cópia
fora do iCloud, porque os diretórios de tempo ocupam cerca de 0,5 GB:

```sh
cp -R case /tmp/marangoni && cd /tmp/marangoni
openfoam ./Allrun                 # malha, recorte do barco, 8 processos, reconstrução (~10 min)
python3 /caminho/export_fields.py /tmp/marangoni -o /caminho/data/marangoni_fields.npz
```

## Modelo

- Janela de 16 × 8 × 2 cm de água (ν = 10⁻⁶ m²/s) que **acompanha o barco** (referencial do barco),
  malha 128 × 64 × 24 refinada perto da superfície. O casco de 3 × 1,5 cm com 3 mm de calado é
  recortado (`subsetMesh`) e fica como parede.
- Movimento: `boatMotion` (function object `coded`, em `system/controlDict`) soma a cada passo
  a força de linha de contato, a de pressão e a viscosa no casco, e integra m dU/dt = Fₓ com
  m = 1,35 g (água deslocada). Depois impõe −U(t) na entrada (proa) e no fundo da cuba e publica
  a = dU/dt, que entra como força fictícia −a na água (`constant/fvOptions`).
- Superfície livre plana (Ca ≈ 10⁻⁴), com condição `slip`.
- Tinta: escalar passivo `c` (`scalarTransport`), D = 2×10⁻⁶ m²/s (efetivo, numérico).
  A fonte presa à popa relaxa c → 1 por 3 s e depois desliga (tinta finita).
- γ = γ₀ − Γc, com Γ = 1 mN/m. A tensão Marangoni τ = ∇ₛγ entra como força de corpo
  f = τ/(ρΔz) na camada de células da superfície.
- O app converte para o referencial do laboratório: x_lab = x + X(t), u_lab = u + U(t).

**Resultado:** com tinta, U chega a 1,1 cm/s (a ordem de cm/s do experimento); depois decai
sob o arrasto. O barco percorre 5,9 cm em 8 s.

**Fora do modelo:** guinada e arfagem; deformação da superfície; surfactante só de superfície
com adsorção; γ(c) não linear. O cisalhamento no casco omite a parte normal do tensor viscoso
(~2% de diferença em relação ao `forces`). Os números são ilustrativos; para um barco mais
rápido, aumente Γ em `constant/fvOptions`, `system/controlDict` e `export_fields.py`.
