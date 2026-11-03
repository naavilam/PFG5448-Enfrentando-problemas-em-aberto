# Ink Boat — Marangoni Propulsion (reservatório finito)

Execute `pip install -r requirements.txt` e `streamlit run app.py` a partir desta pasta.

O modelo mantém transporte superficial de uma espécie efetiva, resposta espectral de Stokes, força capilar na borda circular e arrasto aproximado. Acrescenta um reservatório de tinta finito de volume inicial `V0` (µL), vazão volumétrica `Q` (µL/s), densidade da tinta e fração efetiva `eta` da massa injetada que entra na concentração superficial Γ. A cada subpasso numérico, `released = min(remaining, Q*dt)`, conservando exatamente o orçamento de tinta. A injeção para quando o reservatório esvazia. O campo de concentração já presente permanece e relaxa por transporte, difusão e remoção.

A equação de Γ é uma aproximação de espécie superficial, não um VOF trifásico de tinta imiscível. A taxa Q constante é uma hipótese, não uma lei experimental; beta e eta requerem calibração independente. O esgotamento não força a velocidade a zero imediatamente nem garante um único máximo. As tensões tinta-ar e tinta-água só entram no coeficiente de espalhamento diagnóstico. A resposta de Stokes requer Re pequeno, o que deve ser verificado nos resultados. A solução usa contorno periódico no plano e não resolve a presença do corpo no campo de escoamento.

Mudanças em V0, raio, malha e domínio exigem clicar em Reiniciar.
