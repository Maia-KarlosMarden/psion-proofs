# P14 / WP-BIF: bifurcação de nível na formação (pré-registro)

Hipótese (engenharia reversa de WP-CAL): cada trio formado escolhe o degrau B ou L
(setor dual de §6), com viés ε de sinal oposto para matéria e antimatéria; depois da
formação o degrau é conservado (rótulo de P14-B, aqui POSTULADO). Pergunta: B = −L,
Q = 0 e "hidrogênio + 2ν̄" saem da estrutura sem imposição? Quanto vale η/ε em função
do ritmo da montagem? Que combinação ε·κ·D a PH precisaria entregar?

## Rede (base: WP-NOW etapa A, ontologia O−)
Espécies: 6 ψ, 6 ψ̄, e para cada uma das 20 tuplas t os quatro estados t_B, t_L, t̄_B,
t̄_L (92 espécies). Janela 1,38e11 → 1e8 GeV, Γ/H de WP-NOW, K = 10, pesos L1/L2, g = 1.
- Formação/dissociação por canal: d t_r/dx ∋ ½ k_t P_aP_bP_c − k_t t_r/K (r = B, L);
  consumo de psions por J_t = k_t (P_aP_bP_c − (t_B + t_L)/K).
- Viés de CP (forma de leptogênese, sem assimetria em equilíbrio): matéria
  d t_B ∋ +½ ε s_t J_t, d t_L ∋ −½ ε s_t J_t; antimatéria com −ε (CP).
- Transições t → t′ (W de WP-NOW) e rearranjos 2→2 dentro de cada degrau; o degrau
  é conservado por todo processo depois da formação.
- Aniquilação t_r t̄_r ↔ radiação (equilíbrio (Y_eq,t/2)²); ψψ̄ como em WP-NOW.
- Pesos do viés: S1 s_t = 1 (orientação universal); S2 s_t = (n_L − n_R)/3 (quiralidade
  do trio em relação ao modo zero de J.6).

## Grade
c_Γ ∈ {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1} × lei ∈ {L1, L2} × S ∈ {S1, S2}, ε = 1e-3
(linearidade: ε = 1e-2 no melhor ponto). 24 corridas.

## Saídas
B = Σ(t_B − t̄_B), L = Σ(t_L − t̄_L), Δ_livre, N_ψ, Q total; η = 7,04 Y0 B (Y0 de WP-SIM);
η/ε; ε necessário para η = 6,1e-10 com diluição D = 1 e D = ξ*⁻³ (ξ* = 3,6 e 39, P14-c).
Etapa final (contabilidade, não simulação): com B, L, Q conservados e os estados mais
leves p, e, ν: N_p = B, N_e = N_p (neutralidade), N_ν̄ − N_ν = −L − N_e.

## Controles
C1 ε = 0: B = L = 0 (< 1e-12). C2 |N_ψ| < 1e-12. C3 B + L = 0 (< 1e-12 relativo a |B|).
C4 Q total = 0 (< 1e-12). C5 equilíbrio com taxas constantes (r = 1e3, tempo 20), ε = 0,1:
assimetria < 1e-12. C6 ε → −ε: B → −B (1e-6 relativo).

## Critério e antecipação registrada
Estrutura confirmada se C3, C4 passam e a contabilidade final dá 2 ν̄ líquidos por
hidrogênio. Espero: tudo isso sai exato; η/ε máximo perto da transição (c_Γ ~ 1e-3–1e-2),
de ordem 1e-6–1e-5, e nulo nos dois extremos; S2 menor que S1 por cancelamento entre
trios L-ricos e R-ricos. Isto testa só a estrutura; o número exige derivar ε, κ, D.
