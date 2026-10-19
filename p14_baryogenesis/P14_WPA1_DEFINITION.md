# P14 / WP-A, passo 1: qual transição, qual potencial, qual ξ* é possível (pré-registro)

Entradas proibidas: η_B, ξ* = 434, Ω_DM.

## Achados de leitura (antes de rodar)

1. §13.1(3) e K.9 afirmam "transição Z6 de primeira ordem em T = Λmeta" sem
   potencial, sem calor latente, sem ordem de grandeza do super-resfriamento.
2. J.7 V.5: o sistema estatístico relevante é Z2 × Z3, cada fator com transição
   de Kramers–Wannier autodual (picture tipo Ising × Potts de 3 estados).
3. QFT_PIPELINE.md calibra o gap NJL a T = 0 com σ₀ = m_Z6 = 2,9961 MeV e corte
   Λmeta = 1,26×10¹² GeV, e chama a nucleação em T_nuc = Λ/ξ* de "primeira
   ordem" (condição de Maxwell).
   Tensão: um condensado com gap de 3 MeV tem T_c de ordem MeV, não 10⁹–10¹² GeV.
   A transição da bariogênese e o condensado quiral de 3 MeV não podem ser o
   mesmo fenômeno.

## Fato geral que organiza o teste

ξ* = T_c/T_nuc nunca excede T_c/T_0, onde T_0 é a temperatura espinodal em que a
fase simétrica deixa de ser mínimo local: abaixo de T_0 não há barreira, e a
transição ocorre por rolamento, não por nucleação. Logo cada candidato a
potencial dá um teto para ξ* sem precisar do bounce.

## Modelos candidatos (os três que o texto sugere), todos computados

M1 — Condensado quiral NJL em campo médio, corte Λmeta, N_c = 6 estados
     internos, N_D ∈ {2,4}, G fixado por σ₀ = 2,9961 MeV. Calcular V(σ,T),
     a ordem da transição (sinal do coeficiente quártico de Landau em T_c),
     T_c, T_0 e o teto T_c/T_0.
M2 — Potts de q = 3 estados em campo médio (o fator Z3 de V.5). Analítico:
     T_c = zJ(q−2)/[2(q−1) ln(q−1)], T_0 = zJ/q. Teto T_c/T_0. Conferir
     numericamente minimizando a energia livre.
M3 — Potencial classicamente conforme (Coleman–Weinberg) com massa térmica,
     V = (λ/4) σ⁴ [ln(σ²/Λ²) − 1/2] + c T² σ²/2. Aqui T_0 → 0 e não há teto:
     ξ* é fixado pela ação de bounce. Reportar ξ*(λ) numa grade de λ, só para
     mostrar que faixa de λ daria a faixa-alvo de WP-C (44 a 302). **M3 não está
     no paper**: usá-lo é um postulado novo (P14-c), e o mapeamento de g² = 2π²/3
     e S em λ é escolha a declarar. Neste passo M3 é exploratório e não conta
     como resultado da PH.

## Critério pré-registrado

Se M1 e M2 derem teto ≤ 2, fica registrado que **nenhuma transição descrita no
paper permite ξ* na faixa-alvo**, e que ξ* ~ 10²–10³ exige um potencial quase
conforme que o paper não contém. Isso é um achado negativo sobre a eq. (33)
como está, não sobre a PH em geral. Se M1 ou M2 derem teto > 10, seguir para o
bounce no passo 2.
