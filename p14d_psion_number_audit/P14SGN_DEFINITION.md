# P14 / WP-SGN: a fase do vácuo Z6 é uma fase de CP física no setor de triplas? (pré-registro)

Condição C2 do artigo: o viés da bifurcação precisa de sinal fixado pela orientação do
condensado. Pergunta aqui, calculável com operadores do paper: a fase do vácuo Z6,
w_k = e^{iθ_k}, θ_k = 2πk/6 (K.10; WP-CP A3), posta no termo de enrolamento de J.6.2,
produz uma fase de CP FÍSICA (não removível por refaseamento dos 6 modos) no
Hamiltoniano das 20 triplas? É ímpar sob k → 6−k? Nula em k = 0, 3?

## Modelo (base de Fock de WP-J1 / WP-SIM, 6 modos, 20 triplas)
H_k = h + g O_JS + λ V_k, com h = diag|m| + bond t = 1/2 em |2⟩↔|3⟩ (J.6.3), O_JS de
WP-SIM, g = 1, λ = 1/6 (fração de enrolamento), e
V_k = e^{iθ_k} (J⁺P_L)(J⁺P_R) + e^{−iθ_k} (J⁻P_L)(J⁻P_R)   (J⁺Σ⁺ + J⁻Σ⁻ com fase Z6).
O paper não fixa onde a fase se acopla ao bond gerado pelo mesmo condensado; três
colocações declaradas: PA bond t e^{iθ_k}, V com e^{iθ_k}; PB bond t e^{iθ_k}, V com
e^{2iθ_k}; PC bond real, V com e^{iθ_k}. Antimatéria: H_k*.

## Medidas
M1 Invariante de triângulo: J3(k) = max_{a,b,c} |Im(H_ab H_bc H_ca)| na base de triplas,
   e com sinal: J3s(k) = Im(H_ab H_bc H_ca) no triângulo que maximiza |·| em k = 1.
M2 Removibilidade: min sobre α ∈ R⁶ de ‖Im(U_α H_k U_α†)‖_F, U_α = diag(e^{iα_k}) induzido
   no espaço de Fock (50 partidas aleatórias, Nelder–Mead + BFGS). Física se > 1e-8.
M3 Assimetria dinâmica: rede quântica de WP-LOOP com H_k no lugar da fase inserida,
   ω = 10, γ = 1, c_Γ = 1e-3, lei L2 (melhor ponto de WP-LOOP); saída Δp(k), Δe(k),
   para k = 0..5.

## Controles
C1 k = 0 e k = 3: J3 = 0, M2 = 0, Δ = 0 (< 1e-12). C2 J3s(6−k) = −J3s(k) e Δ(6−k) = −Δ(k)
(1e-6 relativo). C3 H_k hermitiano (< 1e-12).

## Antecipação registrada
PA: fase removível (bond e V coerentes, mesmo número de inserções → α_k = −kθ/2 remove
ambos? não: o bond tem Δm = 1 e V tem Δm = 2; com α_k linear o bond ganha b e V ganha 2b;
PA fica com fase física θ/2; PB (bond θ, V 2θ) é removível; PC física θ). Logo espero:
PA e PC físicas com J3 ≠ 0 e ímpar, PB removível com J3 = 0 para todo k. Em k = 0, 3,
tudo zero. Δp(k) ímpar em k, de ordem 1e-3–1e-2 (fase de ordem 1, contra 1e-3 em WP-LOOP).
Isto fixa o SINAL por vácuo e mostra que a fonte de CP existe na PH (espontânea). O
TAMANHO do viés de bifurcação continua fora do alcance sem C1.
