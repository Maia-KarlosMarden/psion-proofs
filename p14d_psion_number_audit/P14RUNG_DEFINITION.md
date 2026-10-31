# P14 / WP-RUNG: seleção por proto-cor e degrau como enrolamento (pré-registro)

## Parte 1 (estrutural, sem hipótese nova)
Proto-cor do rótulo: c(k) = k mod 3 (J.6.4: Q_k = q(k mod 3) + ⌊k/3⌋). Contar as triplas
de proto-cor completa {0,1,2} (singletos) e verificar se são exatamente os estados
identificados no paper (e, p, n/ν_τ, Σ⁰/ν_μ, Λ⁰/ν_e e as cópias carregadas positivas),
e se as 12 restantes (não singletos) incluem as duas neutras com D = 2 e as 8 de carga
fracionária. Ligação com a "sombra de 8 estados" de App. A.4.

## Parte 2 (hipótese H-W, declarada)
H-W: o degrau ℓ é o número de quanta de enrolamento do condensado presos ao compósito;
cada quantum custa ação 1/b_e (daí m(ℓ) = Λ e^{−ℓ/b_e}) e, no vácuo k, carrega a fase
θ_k = 2πk/6. O degrau é então uma carga de enrolamento (π₁ = Z de WP-TOP), só alterável
passando pelo núcleo desconfinado (C1).
Atribuição de degraus (declarada): B = 8 para toda tripla; L = 10 para triplas carregadas
e L = 15 para neutras (ν: m(15) da escala de Δm²₂₁, §11). Δℓ_t = ℓ_L − ℓ_B ∈ {2, 7}.

## Parte 3 (estimativa do viés, [T4])
Formação em dois canais com amplitudes a_B e^{iℓ_Bθ}, a_L e^{iℓ_Lθ} (|a_B| = |a_L|) e
reespalhamento entre canais com parte absortiva κ: ε_t = κ sin(Δℓ_t θ_k) (forma padrão de
interferência fraca × forte). κ declarado em três valores: α_meta/π = 1/6, g²/16π = 0,131
(g² = 2π²/3), e o enrolamento λ² = 1/36.
Rede de WP-BIF com pesos s_t = sin(Δℓ_t θ_k)/max e ε = κ·max|sin|, só nas 8 triplas
singleto (as não singletos sem viés), lei L1 e L2, vácuos k = 1, 2,
T_RH ∈ {1,38e10, 4e9, 1,4e9, 7e8, 3e8, 1,5e8, 7e7}, janela até T_RH/1e3.
Saídas: η(T_RH); sinal de η por vácuo; T_RH* que dá η = 6,1e-10 (interpolação log).

## Controles
C1 k = 0, 3: η = 0 (< 1e-15). C2 k → 6−k inverte η (1e-6). C3 B + L = 0, Q = 0, N_ψ = 0.

## Antecipação registrada
Parte 1: 8 singletos = conjunto observado; 12 não singletos = exóticos. Parte 3: ε ~ 0,02–0,1,
grande demais no pico; η observado só para T_RH* ~ 1e8–5e8 GeV, abaixo da janela ótima.
Em k = 2, sin(2θ) e sin(7θ) têm sinais opostos (cancelamento parcial entre carregadas e
neutras); em k = 1, mesmo sinal. Tudo isto é condicional a H-W.
