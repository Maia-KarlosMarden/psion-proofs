# P14 / WP-M: massas de tupla e conjugada no mesmo nível de cascata (pré-registro)

Leitura adotada pelo autor (19/10/2026): conjugação C_A (k → 6−k) com o nível de
cascata ℓ como número quântico separado (§8.1 do v85):
  próton = (0,4,5) em ℓ_p, antipróton = (0,1,2) em ℓ_p;
  elétron = (0,1,2) em ℓ_e = 10, pósitron = (0,4,5) em ℓ_e.
Exigência experimental: m(p) = m(p̄) a ~1e-11 (e m(e⁻) = m(e⁺)).

## Três leituras da regra de massa, todas calculadas
- M1: eq. (10) pura, m(ℓ) = Λ exp(−ℓ/b_e). Massa só do nível.
- M2: eq. (10) × eq. (12) com |m_J,R| = |Σ m_J dos constituintes R da tupla|
  (leitura do código `P_yield_ratio.py::_mass`).
- M3: eq. (10) × eq. (12) com |m_J,R| como rótulo de excitação (§6, Tabela 5):
  0 no estado fundamental de cada tupla.
Constantes: Λ = 1,2609e12 GeV, b_e = 325/1152, A = 13/12, p = 5/6.

## Testes
T1. Para ℓ_e = 10 e ℓ_p = 8 (o paper: m(8) ≈ 610 MeV, escala do núcleon),
    massas de (0,1,2) e (0,4,5) em cada leitura; razão m(t)/m(C_A(t)).
T2. Sanidade: em cada leitura, a massa do próton (0,4,5) em ℓ = 8 contra 938 MeV
    (ordem de grandeza; o paper só afirma "escala do núcleon" em ordem dominante)
    e do elétron em ℓ = 10 contra 0,511 MeV.
T3. Contribuição do Hamiltoniano interno de J.6.3 (diag|m_k| + salto t = 1/2):
    ⟨t|h|t⟩ − ⟨C_A(t)|h|C_A(t)⟩ para o par (0,1,2)/(0,4,5), em unidades de h,
    e a fração ímpar de h sob C_A (já medida em WP-J1: 0,26). Registrar que, se h
    contribui para a massa de repouso, a degenerescência exige que essa
    contribuição seja suprimida ou simétrica.

## Critério (pré-registrado)
A leitura passa se m(t) = m(C_A(t)) exatamente no mesmo nível para o par
(0,1,2)/(0,4,5). Reportar quais passam. Se M1 ou M3 passarem, a leitura
"C_A + nível" é compatível com a igualdade de massas, e a assimetria tem de
vir de taxas de formação fora do equilíbrio (próximo passo), não de massas.
