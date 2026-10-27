# Ink Boat — Etapa 1 (protótipo exploratório)

## Execução

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## O que resolve

Modelo de lubrificação para filme fino em substrato rígido, acoplado a uma espécie superficial insolúvel:

- `∂t h + ∇·q = 0`
- `q = -(h³/3μ)∇p + (h²/2μ)∇γ`
- `p = -γ₀∇²h + ρg(h-h₀)`
- `∂t Γ + ∇·(Γu_s) = D_s∇²Γ + S`
- `u_s = -(h²/2μ)∇p + (h/μ)∇γ`
- `γ = γ₀ − EΓ`

**Limitações importantes**

1. **Não representa ainda duas fases móveis tinta–água**: a hipótese de substrato rígido é inadequada para previsões quantitativas do Ink Boat.
2. Não há tensão interfacial tinta–água nem linha de contato móvel explícita.
3. A anisotropia é introduzida na fonte por uma perturbação angular; o modelo não garante crescimento espontâneo de pétalas.
4. Usa diferenças centrais explícitas e clipping de positividade: é um protótipo didático, não um solucionador CFD validado. Para Pe elevado, implementar fluxo upwind/conservativo e testes de convergência.
5. `Γ` representa espécie superficial ativa, **não** diretamente a espessura de uma tinta imiscível. A intensidade ótica não equivale automaticamente a `Γ`.
6. Os parâmetros foram escolhidos para exploração, não ajustados aos experimentos.

## Próxima versão física

Para simular os padrões observados, usar modelo de duas fases com interface móvel e transporte de composição/surfactante, ou uma derivação de filme fino de duas camadas com mobilidade da água. Validar primeiro uma frente radial e depois sua estabilidade angular.
