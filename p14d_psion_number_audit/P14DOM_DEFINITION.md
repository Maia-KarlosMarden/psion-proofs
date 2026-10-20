# P14 / WP-DOM: montagem de compósitos dentro de um domínio Z6 que viola CP (pré-registro)

Pergunta do autor: antes de o compósito existir, a montagem ψψψ → e⁻ pode ser mais
rápida que ψ̄ψ̄ψ̄ → e⁺? Na eq. (1) a CP de campo é exata (WP-CP); dentro de um
domínio de vácuo Z6 ela está quebrada espontaneamente. Testar o que a montagem
produz nesse domínio, respeitando CPT e unitariedade, e o que sobra globalmente.

## Parte A: rede de montagem num domínio (ontologia O−, rede de WP-NOW etapa A)
Mesma rede de 52 espécies, mesma janela 1,38e11 → 1e8 GeV, mesmas taxas de contato,
K = 10, δ = 0 (psions = antipsions no início), pesos de ligação L1/L2, g = 1, φ = 0.
Ligação ψψψ ↔ t idêntica nos dois setores: um trio de psions só forma uma tupla, e
CPT + unitariedade proíbem assimetria num canal único (Γ(ψψψ→t) = Γ(t̄→ψ̄ψ̄ψ̄)
e a constante de equilíbrio é fixada pelas massas iguais).

Violação de CP do domínio, de magnitude ε (a PH não dá o valor):
- Transições t → t′ (arestas do grafo de h + J·Σ): matéria R_m(i→j) = W_ij(1 + ε x_ij),
  x antissimétrico; antimatéria R_a = R_mᵀ (CPT: Γ(t→t′) = Γ(t̄′→t̄)).
- Rearranjos ab ↔ cd: matéria com taxa direta (1 + ε c_r) e inversa (1 − ε c_r);
  antimatéria com os sinais trocados (CPT).
- Unitariedade: o vetor (W x, c) está no núcleo da matriz de incidência (cada tupla
  recebe tanto quanto perde no equilíbrio uniforme), de modo que o equilíbrio
  continua simétrico e qualquer assimetria só nasce fora do equilíbrio.
- Direção do vetor: a PH não a fixa; sorteio gaussiano projetado no núcleo, 5
  sementes (0..4), normalizado a max |x|, |c| = 1.
Variantes: V1 só transições; V2 transições + rearranjos (núcleo conjunto).
Grade: V ∈ {V1, V2} × c_Γ ∈ {1e-5, 1e-3, 1} × lei ∈ {L1, L2} × ε ∈ {1e-6, 1e-3} ×
semente ∈ {0..4}: 120 corridas.
Saídas: Δ_p (tupla (0,4,5) menos antitupla), Δ_e ((0,1,2)), Δ_livre = Σ(ψ_k − ψ̄_k),
N_ψ líquido, e se há "hidrogênio do domínio" (Δ_p > 0 e Δ_e = Δ_p em sinal e ordem).

## Parte B: o domínio existe na época da montagem? E o que sobra globalmente?
- K.10: as paredes Z6 aparecem em T ≈ 3 MeV e se aniquilam em ~4 anos. Calcular
  Γ/H das taxas de montagem em 3 MeV (c_Γ = 1) e se a montagem ainda opera.
- Vácuos conjugados são degenerados (ΔC simétrico): volumes iguais de ε e −ε; a
  média global da assimetria é zero.
- Tamanho comóvel máximo de um domínio hoje = horizonte em T_dom, para T_dom = 3 MeV
  e para a hipótese T_dom = 1,38e11 GeV. Comparar com o limite observacional para
  domínios de antimatéria (≳ Gpc; Cohen, De Rújula e Glashow 1998).

## Controles
C1 ε = 0: assimetrias < 1e-12 (o piso de arredondamento medido em WP-NOW é 2e-13).
C2 N_ψ líquido conservado: |N_ψ| < 1e-12 em toda corrida.
C3 domínio conjugado: Δ(−ε) = −Δ(ε) (desvio relativo < 1e-6).
C4 equilíbrio: com taxas constantes e grandes (r = 1e3, tempo 20) a partir do
   equilíbrio uniforme, a assimetria fica < 1e-12 (a construção não cria assimetria
   em equilíbrio).

## Antecipação registrada
O grafo de transições tem 5 componentes; a do elétron, {(0,1,2), (0,1,3)}, é uma
aresta só, sem ciclo, e não comporta circulação: em V1 espero Δ_e = 0 exato e Δ_p ≠ 0
(o próton está na componente grande, com ciclos). Espero Δ ∝ ε, Δ_livre = −3 ΣΔ_tuplas
(o reservatório de antipsions), média global zero, e que em 3 MeV a montagem esteja
congelada (Γ/H ≪ 1), de modo que domínios de K.10 chegam tarde demais.
