# Ink Boat — implementação reduzida inspirada no Diário de Bordo

Execute `streamlit run app.py`.

## Correspondência com o diário
- Kang (2021), capítulo 3: força e torque sobre o corpo (equações 3.5 e 3.6): **implementados como EDOs 2D**, com força capilar integrada e arrastos efetivos.
- Lei de tensão superficial linear em concentração: **implementada** na borda, usando concentração angular fenomenológica adimensional.
- Navier–Stokes 3D, cisalhamento tangencial de Marangoni, campo de pressão e integral de tensão viscosa na superfície molhada: **não implementados**.
- Transporte 3D de espécie solúvel do artigo: **não implementado**; tinta insolúvel exige modelo distinto.
- Três interfaces: tensões parametrizadas e coeficiente de espalhamento S; **não há interface multifásica móvel resolvida**.
- Re/Pe são diagnósticos estimados, não condições que alteram o solucionador.
- Rotação: a tração radial em um barco circular não produz torque. O controle de braço lateral aplica uma hipótese de força excêntrica para explorar a resposta rotacional.

## Limites e validação
Esta implementação **não é uma previsão ab initio** de velocidade: o contraste de tensão, a fonte angular e os arrastos devem ser medidos ou calibrados. Para um modelo hidrodinâmico fiel, resolver campo de velocidades 3D com a condição tangencial de Marangoni e acoplamento ao movimento do barco; para tinta imiscível, incorporar a interface tinta–água e a linha de contato móvel.
