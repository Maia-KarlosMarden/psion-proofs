# P14 / WP-TPROG: programa de temperatura livre (pré-registro)

Pergunta do autor: diferentes temperaturas favorecem diferentes compósitos.
Teste: com a PH como está (ontologia O−, N_ψ e CP exatos), uma história de
temperatura qualquer (resfriamentos, reaquecimentos, repetições) muda a
composição por nível? Cria assimetria?

## Modelo (setores partícula/antipartícula, um por nível)
Massas da cascata, eq. (10), m(ℓ) = Λ e^{−ℓ/b_e}, com 1/b_e = ln(610/0,511)/2
= 3,5424 e Λ = 1,2387e12 GeV (âncoras ℓ = 8 → 610 MeV, ℓ = 10 → 0,511 MeV):
- S_top: tupla do topo (nível 0,8; m = 7,28e10 GeV; g = 2). ⟨σv⟩ = 1/m²
  (acoplamento forte; entrada declarada, sem valor na PH).
- S8: núcleon, nível 8, m = 0,938 GeV (sensibilidade: 0,610), g = 2,
  ⟨σv⟩ = 1/m_π² (taxa hadrônica medida, a mesma de WP-NOW).
- S10: elétron, nível 10, m = 0,511 MeV, g = 2, ⟨σv⟩ = π α²/m² (e⁺e⁻ → γγ).
O nível 9 (X17) é bóson ψψ̄, autoconjugado (N_ψ = 0): não tem assimetria;
fica fora.

Por setor: Y_X, Y_X̄ (Y = n/s). Processo único: X X̄ ↔ radiação,
dY_X̄/d ln a = −(s⟨σv⟩/H) [Y_X Y_X̄ − Y_eq²], idem Y_X; A = Y_X − Y_X̄
fica constante fora dos reaquecimentos. Y_eq Maxwell–Boltzmann exato,
n_eq = g T³ (x²/2) K₂(x)/π², x = m/T. g*(T) do Modelo Padrão em degraus
(106,75 / 86,25 / 75,75 / 61,75 / 17,25 / 10,75 / 3,91 com cortes em 80, 4,
1, 0,15, 0,1 e 5e-4 GeV). H = 1,66 √g* T²/M_Pl, M_Pl = 1,2209e19 GeV.
Passo implícito em ln a (Euler para trás, solução exata da quadrática),
Δln T = 0,004 (convergência: 0,001).

Reaquecimento de T₁ a T₂ > T₁: injeção instantânea de entropia; todo Y é
multiplicado por s(T₁)/s(T₂). (Uma subida adiabática exigiria contração; fica
fora.)

Saídas: n/n_γ hoje = 7,04 Y (entropia conservada depois do fim), para X e X̄
de cada setor; Ω h² = 2,74e8 (m/GeV)(Y_X + Y_X̄) para S_top (diagnóstico).

## Histórias declaradas (todas terminam em T = 1e-5 GeV)
H1 resfriamento monotônico de 1e12 GeV.
H2 H1 com reaquecimento a 2 GeV quando T = 0,01 GeV.
H3 H1 com reaquecimento a 1e11 GeV quando T = 1e6 GeV (acima do congelamento
   do topo).
H4 "serrote": 1e12 → 0,05 → 5 → 0,05 → 5 → 1e-5 GeV.
H5 H1 com reaquecimento a 0,05 GeV quando T = 1e-4 GeV (acima do congelamento
   do elétron, abaixo do núcleon).
Cada uma com η_in = 0 e com η_in = 1e-9 (posto à mão em S8 e S10, ⟨A⟩ = η/7,04,
neutralidade de carga); S_top sem assimetria.

## Controles
C1 η_in = 0: A = 0 em todo setor e toda história (exato por construção:
   cada aniquilação remove um X e um X̄. Não é teste novo da PH; a rede
   completa já foi testada em WP-NOW).
C2 η_in ≠ 0: A final = A_in × ∏ fatores de diluição (a assimetria nunca cresce).
C3 convergência no passo: Δ ln T 0,004 vs 0,001, diferença < 1%.
C4 implementação JavaScript da página contra este Python: < 1e-6 relativo
   (mesmo algoritmo e passo).

## Antecipação registrada
Assimetria nunca aumenta; reaquecimento acima de m/25 de um setor reapaga a
relíquia simétrica daquele setor e a recongela no mesmo valor (ela depende só
de m e ⟨σv⟩, não da história anterior), enquanto setores mais pesados que
T₂ ficam só diluídos. A composição relativa entre níveis muda com a história;
a razão matéria/antimatéria, não.
