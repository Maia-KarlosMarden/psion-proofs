# P14 / WP-SIM: rede de Boltzmann de formação de tuplas (pré-registro)

Leitura: conjugação C_A (k → 6−k) com nível de cascata separado (WP-M):
p = (0,4,5)@8, p̄ = (0,1,2)@8, e⁻ = (0,1,2)@10, e⁺ = (0,4,5)@10.
Massas iguais no mesmo nível (eq. 10); janela T ≫ massas (relativístico).
Entradas proibidas: η_B (só na comparação final), ξ* = 434.

## Espécies
6 psions livres (Y_k) + 20 tuplas × 2 níveis (Y_{t,ℓ}), ℓ ∈ {8, 10}: 46 espécies.
Y em unidades de n_ψ/s inicial (Σ Y_k = 1 no início).

## Processos (todos com balanço detalhado; massas iguais ⇒ equilíbrio simétrico)
P1 ligação/dissociação: ψ_a ψ_b ψ_c ↔ t@ℓ, t = {a,b,c}; fluxo
   J = k_b(t) u_ℓ [Y_a Y_b Y_c − Y_{t,ℓ}/K], u_8 = u_10 = 1/2, K = 10.
P2 transições internas 1→1 no mesmo nível: t ↔ t′ com peso W = |⟨t′|H|t⟩|²,
   H = h + g O_JS (h de J.6.3 com o salto t = 1/2; O_JS na convenção D1, como em
   WP-J1). Inclui o salto 2↔3 do condensado e os deslocamentos J·Σ.
P3 rearranjo 2→2 no mesmo nível (multiconjunto de rótulos preservado), peso 1.
Sem aniquilação (pares relativísticos na janela).

## Lei de taxa (pesos de ligação k_b(t), média normalizada a 1)
L1 regra de ouro: k_b(t) ∝ |⟨t|H|t⟩|².
L2 Boltzmann interno: k_b(t) ∝ exp(−⟨t|H|t⟩) (H em unidades de h).

## Cosmologia
Início em T_RH = 1,38e11 GeV (linha λ = 1/25, c = 1 de WP-A), só psions livres.
Fim em 1e8 GeV. x = ln(T_RH/T); dY/dx = (Γ/H) F(Y), Γ/H = c_Γ g⁴ T³ M_Pl /
(1,66 √g* Λ⁴), g² = 2π²/3, g* = 106,75, Λ = 1,2609e12 GeV.

## Grade (12 corridas)
c_Γ ∈ {1e-5, 1e-3, 1} × g ∈ {0,1; 1} × lei ∈ {L1, L2}.

## Saídas
Δ_p = Y(0,4,5)@8 − Y(0,1,2)@8 (excesso de prótons sobre antiprótons);
Δ_e = Y(0,1,2)@10 − Y(0,4,5)@10 (excesso de elétrons sobre pósitrons);
carga líquida dos compostos; η = 7,04 · y₀ · Δ_p, y₀ = 0,208 g_ψ/g* com g_ψ = 12.

## Controles (devem passar, senão nada conta)
C1 pesos simetrizados sob C_A ⇒ |Δ| < 1e-12.
C2 conservação de N_ψ = Σ Y_k + 3 Σ Y_t: desvio relativo < 1e-8.
C3 identidade dos níveis: |Y_{t,8} − Y_{t,10}| < 1e-12.

## Critério
- Matéria de hidrogênio: Δ_p > 0 e Δ_e > 0 em alguma corrida.
- Magnitude: η a um fator 3 de 6,1e-10 conta só se o sinal for o de hidrogênio.

## Transparência (antecipação ao projetar)
Com u_8 = u_10 e todos os processos cegos ao nível, os dois níveis evoluem como
cópias idênticas (C3), logo Δ_e = −Δ_p exatamente: taxas que dependem só dos
rótulos favorecem a mesma tupla nos dois níveis, e produzem prótons com pósitrons
ou antiprótons com elétrons, nunca hidrogênio. O simulador verifica isso e mede
|Δ_p| e a lavagem em função de c_Γ.
