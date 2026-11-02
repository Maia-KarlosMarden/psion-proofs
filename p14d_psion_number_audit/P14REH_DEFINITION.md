# P14 / WP-REH: reaquecimento gradual, degraus consistentes e normalização (pré-registro)

Responde aos testes 1–3 do parecer de referee (M1–M3) e testa se a diluição do
reaquecimento remove a superprodução (M4).

## Rede
A de WP-COMP (8 singletos × 2 degraus × matéria/antimatéria; criação de pares da radiação
em quatro canais; viés δ = ½ ε s_t ΣJ; E_b = ∞), com:
- Degraus consistentes com WP-UNW (M1): Δℓ = +2 para (0,1,2), +7 para as neutras, −1 para as
  carregadas positivas. s_t = sin(Δℓ_t θ_1)/0,866, ε = κ·0,866, κ = g²/16π = 0,1309, vácuo k = 1.
- Reaquecimento (M2): ρ_φ e ρ_R com dρ_φ/dN = −3ρ_φ − (Γ_φ/H)ρ_φ, dρ_R/dN = −4ρ_R + (Γ_φ/H)ρ_φ,
  N = ln a, H² = (8π/3M_Pl²)(ρ_φ + ρ_R), ρ_R = (π²/30) g* T⁴, g* = 106,75. Γ_φ = H_rad(T_RH).
  Início em N = 0 com ρ_φ = 3H_I²M_Pl²/(8π), ρ_R = 1e-30 ρ_φ, compósitos nulos.
  H_I/Γ_φ ∈ {1e2, 1e4, 1e6, 1e8} (fixa T_max/T_RH; reportar T_max).
- Cinética: dY/dN = (Γ/H) g(Y) − Y d ln(s a³)/dN, com Γ = (T/T_eq)³ H_rad(T) c_Γ e g(Y) os
  termos de colisão de WP-COMP com taxa unitária. Fim quando ρ_φ/ρ_R < 1e-8 e T < T_RH/1e3.
- Normalização (M3): c_Γ ∈ {0,1, 1, 10}.
Grade: T_RH ∈ {1e8, 3e8, 1e9, 3e9, 1e10, 2e10, 4e10, 1e11} × H_I/Γ_φ (4) × L1, L2 × c_Γ (3).
192 corridas. η = 7,04 B.

## Saídas
η(T_RH) por (H_I/Γ_φ, lei, c_Γ); todas as raízes |η| = 6,1e-10 com o sinal; T_max; fração
da grade com |η| entre 1e-10 e 1e-9 (sem sintonia, isso deveria ser uma faixa larga).

## Controles
C1 ε = 0: |η| < 1e-15 (limiar relativo ao sinal típico 1e-6). C2 nos pontos com |B| > 1e-12:
|B + L|/|B| < 1e-9 e |Q|/|B| < 1e-9. C3 conservação da energia no setor ínflaton + radiação:
desvio de (ρ_φ + ρ_R) a³ᐟ⁴-escala esperado; checar que T_RH efetivo (H = Γ_φ) está a ±30% do
nominal.

## Antecipação registrada
Com reaquecimento gradual a produção ocorre perto de T_max e é diluída pela entropia até
T_RH. Espero η reduzido em relação a WP-COMP por fatores grandes quando T_max ≫ T_RH, com
uma faixa mais larga de T_RH dando |η| ~ 1e-10–1e-9 (M4 aliviado). Espero sinal ainda
dependente de lei e T_RH (M1). T_RH* desloca por ~c_Γ^{−1/3} (M3).
