# P14 / WP-TRH: bifurcação de nível em função da temperatura de reaquecimento (pré-registro)

Declaração: antes deste registro rodei UMA checagem exploratória (P14BIF_TRH_exploratory.txt:
c_Γ = 1, ε = 1e-3, S1, L1, cinco valores de T_RH). Este teste a formaliza e a estende.

Pergunta: com acoplamento natural (c_Γ = 1, isto é, G_UV ~ g²/Λ²), o desequilíbrio da
formação é fixado pela temperatura de reaquecimento T_RH? Qual a eficiência η/ε em função
de T_RH e qual ε reproduz η = 6,1e-10 sem diluição adicional?

## Rede
A de WP-BIF (92 espécies, viés S1, ε = 1e-3), com a janela de formação de T_RH até
T_END = 1e8 GeV e Γ/H = c_Γ G4 T³ M_Pl/(1,66 √g* Λ⁴), c_Γ = 1, Λ = 1,2609e12 GeV.
Grade: T_RH ∈ {1,38e11, 7e10, 4e10, 2e10, 1,38e10, 7e9, 4e9, 2e9, 1,4e9, 7e8} GeV ×
lei ∈ {L1, L2}. 20 corridas.

## Saídas
η/ε e ε necessário = 6,1e-10/(η/ε) por ponto; T_eq (Γ/H = 1 com c_Γ = 1); faixa de T_RH
onde ε necessário ≤ 1e-5; T_RH do máximo.
Equivalência: (c_Γ = 1e-3, T_RH = 1,38e11) contra (c_Γ = 1, T_RH = 1,38e10), diferença relativa.

## Controles
C1 |B + L|/|B| < 1e-9 onde |B| > 1e-12. C2 |Q| < 1e-12. C3 |N_ψ| < 1e-12.

## Antecipação registrada
Equivalência a < 1%. Máximo de η/ε ≈ 3,6e-3 perto de T_RH ≈ 1,4e10 GeV; ε necessário
≤ 1e-5 para T_RH ≈ 2e9–4e10 GeV; lavagem acima de ~7e10 GeV. Toda a janela útil fica
abaixo de Λ_meta.
