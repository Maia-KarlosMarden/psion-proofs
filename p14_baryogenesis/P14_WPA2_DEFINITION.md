# P14 / WP-A, passo 2: bounce no potencial de Coleman–Weinberg (pré-registro)

Postulado P14-c (adotado por decisão do autor em 19/10/2026): o parâmetro de
ordem da transição Z6 em Λmeta tem potencial classicamente conforme

    V(φ,T) = (λ/4) φ⁴ [ln(φ²/v²) − 1/2] + (c/2) T² φ²,   v = Λmeta.

Mínimo a T = 0 em φ = v, energia de vácuo da fase simétrica ΔV = λ v⁴/8.
c é o coeficiente de massa térmica; valor de referência c = 1, sensibilidade
com c = 0,3 e c = 3. Entradas proibidas: η_B, ξ* = 434, Ω_DM.

## Candidatos a λ (propostas, não derivações; fixados antes de rodar)

    λ₁ = 1/(2j)²  = 1/25   = 0,0400   (o mesmo 1/(2j)² de ε_ψ)
    λ₂ = ε_ψ      = 25/169 = 0,1479
    λ₃ = 1/N_c    = 1/6    = 0,1667
    λ₄ = b_e      = 325/1152 = 0,2821
    λ₅ = N_c g²/(4π)² = 6·(2π²/3)/(16π²) = 1/4 = 0,2500

Nenhum é derivado do Lagrangiano; são as combinações adimensionais simples da
PH. Reportar ξ* para todos. O λ que reproduzir a faixa-alvo de WP-C (44–302)
é reportado como AJUSTE, não como previsão.

## Cálculo

1. T_c: V(φ_min(T), T) = 0.  Espinodal T_0 = 0 (conforme).
2. Bounce O(3): φ'' + (2/r)φ' = ∂V/∂φ, overshoot/undershoot; S_3(T) =
   4π ∫ r² [φ'²/2 + V(φ) − V(0)] dr.
3. Nucleação: Γ/H⁴ = 1, Γ = T⁴ (S_3/2πT)^{3/2} e^{−S_3/T},
   H² = [ (π²/30) g* T⁴ + ΔV ] / (3 M_P²), M_P reduzida = 2,435e18 GeV.
   Registrar T_eq onde ΔV = ρ_rad (início da dominação do vácuo). Se
   T_nuc < T_eq há inflação térmica; a assimetria produzida é diluída pelo
   reaquecimento até T_RH ≈ (30 ΔV/(π² g*))^{1/4}, e a "ξ*" relevante para a
   eq. (33) deixa de ser T_c/T_nuc. Reportar T_nuc, T_eq, T_RH e ξ* = T_c/T_nuc.
4. Critério: compatível se 44 ≤ ξ* ≤ 302 E T_nuc > T_eq (sem inflação térmica).

## Contagem de tentativas
WP-A: 5 valores de λ × 3 valores de c = 15 tentativas, todas reportadas.
