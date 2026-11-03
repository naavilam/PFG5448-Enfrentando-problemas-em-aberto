# Ink Boat — nova implementação do balanço de forças

Interface Streamlit com equações, constantes e parâmetros à esquerda e gráficos
à direita. Esta implementação é independente de `inkboat_propulsion_finite_ink`;
o código anterior e o PDF da tese são preservados.

## Executar

Na pasta `inkboat_new_implementation`:

```sh
python3 -m pip install -r requirements.txt
python3 -m streamlit run app.py
```

Se usar o venv do projeto, execute esses comandos com o Python desse ambiente.
Nenhum serviço precisa ser publicado. No topo da barra lateral:

- **Iniciar / Retomar:** começa em t=0 ou continua do instante pausado.
- **Pausar:** congela o tempo e mantém as curvas.
- **Parar / Reiniciar:** retorna a t=0, limpando o histórico exibido.
- **Velocidade de reprodução:** controla quantos segundos simulados avançam por segundo real.

A execução para ao atingir a janela de tempo; Iniciar novamente começa outro
ensaio. Alterar parâmetros físicos reinicia o ensaio para não misturar modelos.
Os gráficos e indicadores avançam com o tempo simulado, avaliando a solução
analítica até o instante atual. O download exporta somente o trecho já exibido.
O contraste Δγ(t) gera a força propulsiva; com Δγ=0, o relógio avança mas o barco permanece parado.

## Equações da tese

Saeed Jafari Kang, *Marangoni Propulsion of Active Particles* (2021), §3.2,
páginas numeradas 26–27, PDF nesta pasta.

- Eq. (3.5): ∮ γ t dℓ + ∫ n·σ dS = m dU/dt, U(0)=0.
- Eq. (3.3): σ = −pI + μ[∇u + (∇u)ᵀ].

**ℓₚ é a linha de contato; o parâmetro de comprimento é P = ∮dℓ.**
A direção t da tese é tangente à interface água–ar e normal à linha de contato.
Ela não é a tangente ao próprio contorno; para o disco é a direção radial.
Perímetro sozinho não determina a força: é preciso a distribuição de γ e as
direções ao longo do contorno. Para um contorno fechado plano e γ uniforme,
∮ γ t dℓ = 0. γP é a soma dos módulos, não a força resultante.

## Fechamento aceito para esta etapa

- Disco fino horizontal; translação apenas em x; sem rotação e com imersão fixa.
- P ajustável, R=P/(2π), área inferior molhada A=πR²; borda molhada desprezada.
- Água pura a 20 °C, aproximadamente 0,1 MPa; constantes μ=1,0016×10⁻³ Pa·s
  e ρ≈998,2 kg/m³. γ₀ é ajustável e prescrito, não uma constante inferida.
- Barco e água inicialmente em repouso. γ varia no espaço e no tempo: γ(θ,t)=γ₀+[Δγ(t)/2]cosθ.
- Δγ=γfrente−γtrás; frente é +x, traseira é −x. O contraste satisfaz |Δγ|≤2γ₀ para manter γ≥0. O contraste inicial Δγ₀ é prescrito; permanece constante enquanto há tinta e zera abruptamente no esgotamento, sem resolver o transporte.
- Aproximação de cisalhamento: uₓ(z,t)≈U(t)(1+z/h), para −h≤z≤0, com água em
  repouso em z=−h e não deslizamento na face inferior z=0.
- h é uma escala de cisalhamento efetiva ajustável, não a profundidade física
  conhecida do tanque ou uma espessura calculada da camada limite.
- Pressão sem gradiente horizontal; nenhuma força horizontal de pressão no
  disco simétrico. A contribuição viscosa é Fh≈−μAU/h=−bU.

É um fechamento reduzido assumido, não a solução de Navier–Stokes da tese.
μ não determina sozinha σ. Não se calculam campos de pressão, escoamento
tridimensional, massa adicionada, memória hidrodinâmica nem transporte de tinta.
A aproximação quase estacionária requer atenção aos efeitos de borda, à escala h/R
e ao tempo de difusão viscosa h²/ν em relação ao tempo de resposta m/b. A aplicação
exibe esses diagnósticos e Re; não afirma validade quantitativa em todo o intervalo
dos controles. Peso e empuxo pertencem ao equilíbrio vertical, aqui prescrito.

### Constantes e fonte

[IAPWS SR6-08(2011), propriedades da água líquida a 0,1 MPa](https://iapws.org/technical-guidance/release/LiquidWater),
seção 3, Eq. (7). A correlação de viscosidade é

μ/(10⁻⁶ Pa·s) = 280,68(T/300 K)⁻¹·⁹ + 511,45(T/300 K)⁻⁷·⁷
+ 61,131(T/300 K)⁻¹⁹·⁶ + 0,45903(T/300 K)⁻⁴⁰.

A 293,15 K, fornece aproximadamente 1,0016 mPa·s. Fixamos esse valor para esta
primeira etapa; a temperatura não é um controle variável.

## Gota finita e corte abrupto do contraste

Os controles são volume inicial V₀ (µL), taxa constante de ejeção Q (µL/s)
e contraste inicial Δγ₀. Internamente, 1 µL = 10⁻⁹ m³.

- V(t)=max(V₀−Qt,0).
- t_fim=V₀/Q, para Q>0.
- Δγ(t)=Δγ₀ enquanto V(t)>0; Δγ(t)=0 quando V(t)=0.
- γ(θ,t)=γ₀+[Δγ(t)/2]cosθ.
- Fγ,x(t)=PΔγ(t)/4 e Fγ,y=0.

O volume diminui linearmente, mas o contraste e a força ficam constantes até
t_fim, quando caem abruptamente para zero. Maior Q antecipa esse corte; maior
V₀ o adia. A tensão média γ₀ fica fixa. Q=0 mantém volume e contraste; V₀=0
significa contraste/força zero desde o início. O corte é uma hipótese prescrita,
sem transporte/difusão ou ação residual da tinta já liberada.

A resposta de m dU/dt=Fγ(t)−bU é analítica por trechos:

- Antes de t_fim: U(t)=(F₀/b)(1−exp(−t/τ)), τ=m/b, F₀=PΔγ₀/4.
- Após o corte: U(t)=U(t_fim)exp[−(t−t_fim)/τ].
- A aceleração é [Fγ(t)−bU(t)]/m e salta em −F₀/m no corte.

Velocidade e posição permanecem contínuas: o barco não para instantaneamente.
Os gráficos mostram o volume decrescendo e o contraste em degrau, marcando o
esgotamento com linha pontilhada. A malha de visualização inclui o instante do
corte e seu limite à esquerda para não desenhar uma rampa artificial na força.
O CSV inclui volume, contraste e demais grandezas em SI.

Os valores padrão V₀=0,02 µL e Q=0,005 µL/s dão t_fim=4 s. Δγ₀=0,01 mN/m
é ilustrativo, não medido. Não há força externa.

## O que falta para modelar a experiência

Forma e dimensões reais, massa, imersão/área molhada, temperatura, paredes e fundo
do recipiente, condições iniciais do fluido e uma lei γ(s,t). Para usar diretamente
a Eq. (3.3), é necessário resolver os campos p e u com condições de contorno.

A função `gamma_on_contour(theta, time_s, p)` mantém o contraste inicial
enquanto há tinta e o desliga no esgotamento. Uma futura lei de concentração poderá substituir essa
relação empírica. Um gradiente de tensão real também gera escoamento de Marangoni
na água; esse escoamento não é calculado pelo arrasto −bU e exigirá rever o
fechamento hidrodinâmico. Esta etapa calcula a força capilar direta com contraste
prescrito pelo volume, não a solução completa da tese.

## Verificação

```sh
python3 -m unittest -v test_model.py
```

Os testes cobrem cancelamento para γ uniforme, integral com harmônico conhecido,
repouso sem força, balanço de forças/velocidade terminal, força negativa, limite
sem arrasto, estabilidade em tempos pequenos e conversões geométricas.

Testes adicionais verificam conservação de volume, tempo de esgotamento, contraste em
degrau, continuidade de velocidade/posição, salto da aceleração e relaxação após esgotamento e concordância com uma integração
RK4 independente.

## Streamlit Community Cloud

A pasta atual é `16-ink-boat_WORKING/inkboat_propulsion_finite_ink`, renomeada de
`inkboat_new_implementation`. O entrypoint é
`16-ink-boat_WORKING/inkboat_propulsion_finite_ink/app.py`.

O ambiente de compatibilidade usa Streamlit **1.56.0**, anterior à migração do
servidor Tornado para Starlette/Uvicorn na versão 1.57. O modelo e a interface
foram testados nessa versão em um venv isolado. NumPy 2.4.6, Matplotlib 3.10.9 e
PyArrow 24.0.0 permanecem fixos. O downgrade é um teste de compatibilidade;
a causa dos erros HTTP 500/404 na nuvem não foi comprovada localmente.

Depois de enviar as alterações ao GitHub, reinicie o app em Manage app → Reboot
app. Confirme **streamlit==1.56.0** no log. Ao acessar a aplicação, a linha
`INKBOAT_BOOT` registra a versão efetivamente carregada na sessão Python.

Existe uma rota de diagnóstico sem modelo, NumPy, Matplotlib ou atualização
periódica: abra `https://inkboat.streamlit.app/?diagnostic=1`.

- Se abrir, clique em “Testar comunicação com o servidor”: o contador confirma
  comunicação entre navegador e Python.
- Se continuar apenas carregando e não houver `INKBOAT_BOOT` no log, a sessão
  não chegou ao script. Investigue servidor, frontend estático, proxy ou a
  plataforma antes de alterar as equações.
- Se houver `INKBOAT_BOOT`, mas a página normal falhar, consulte o traceback
  posterior; o diagnóstico separa o carregamento mínimo das importações e gráficos.

A rota não publica segredos nem muda a simulação. Apague `?diagnostic=1` para voltar.
