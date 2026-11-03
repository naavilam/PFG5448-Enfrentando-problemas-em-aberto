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
- Δγ=γfrente−γtrás; frente é +x, traseira é −x. O contraste satisfaz |Δγ|≤2γ₀ para manter γ≥0. O contraste inicial Δγ₀ é prescrito; sua redução é proporcional ao volume restante, sem resolver o transporte.
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

## Gota finita e contraste decrescente

Dois controles definem a duração da propulsão: volume inicial V₀ (µL) e taxa
constante de ejeção Q (µL/s). Internamente são convertidos para m³ e m³/s
(1 µL = 10⁻⁹ m³).

- V(t)=max(V₀−Qt,0).
- t_fim=V₀/Q, para Q>0.
- Δγ(t)=Δγ₀ V(t)/V₀, para V₀>0.
- γ(θ,t)=γ₀+[Δγ(t)/2]cosθ.
- Fγ,x(t)=PΔγ(t)/4 e Fγ,y=0.

O contraste começa no valor inicial e decresce linearmente até zero. Maior Q
esgota a tinta mais cedo; maior V₀ prolonga a propulsão. O valor médio γ₀ não vai
a zero. Q=0 mantém a tinta e o contraste inicial; V₀=0 implica contraste e força
nulos desde t=0. Esses são limites explícitos desta lei fenomenológica.

A hipótese de proporcionalidade ao volume restante é uma escolha deste modelo,
não uma lei derivada de transporte, difusão ou adsorção. A tinta já liberada não
mantém contraste após o reservatório esvaziar nesta aproximação. Δγ₀=0,01 mN/m,
V₀=0,02 µL e Q=0,005 µL/s são valores ilustrativos: t_fim=4 s.

O movimento obedece m dU/dt=Fγ(t)−bU, b=μA/h e U(0)=0. A resposta é calculada
analiticamente por trechos, sem erro de integração temporal: força em rampa
até t_fim e apenas relaxação viscosa depois. O gráfico inclui exatamente o instante
de esgotamento, mesmo se ele cair entre os pontos da malha de visualização.

Após o esgotamento, U(t)=U(t_fim) exp[−(t−t_fim)/τ], com τ=m/b. O barco não para
instantaneamente: ele conserva velocidade e desacelera por arrasto. A velocidade
e a aceleração não precisam diminuir juntas desde o início; U parte de zero,
cresce e depois cai, enquanto a aceleração pode mudar de sinal antes de t_fim.
Contraste negativo inverte o sentido da trajetória.

A interface mostra volume restante/ejetado, contraste atual, forças, velocidade
e aceleração. A linha vertical pontilhada marca t_fim. O CSV inclui volume,
contraste e demais grandezas em unidades SI. Não há força externa.

## O que falta para modelar a experiência

Forma e dimensões reais, massa, imersão/área molhada, temperatura, paredes e fundo
do recipiente, condições iniciais do fluido e uma lei γ(s,t). Para usar diretamente
a Eq. (3.3), é necessário resolver os campos p e u com condições de contorno.

A função `gamma_on_contour(theta, time_s, p)` usa o contraste calculado pela
fração de tinta restante. Uma futura lei de concentração poderá substituir essa
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

Testes adicionais verificam conservação de volume, tempo de esgotamento, contraste
linear, continuidade, relaxação após esgotamento e concordância com uma integração
RK4 independente.
