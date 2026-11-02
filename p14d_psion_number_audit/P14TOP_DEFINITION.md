# P14 / WP-TOP: o número bariônico pode ser topológico ou de simetria na eq. (1)? (pré-registro)

WP-LEV1 mostrou que a regra de níveis precisa conservar B (classe de nível 8)
separado de L, a 1e-8 (contato em Λ_meta) e a 5e-30 (saltos na mesma tupla).
Candidatos clássicos para um B conservado: (i) carga topológica pontual de um
condensado (skyrmion, π₃ da variedade de vácuo, como na QCD); (ii) carga de Noether
contínua; (iii) carga discreta; (iv) rótulo interno do estado.

Declaração: antes deste registro fiz UM cálculo exploratório (simetria contínua da
interação J·Σ nos 6 rótulos: dimensão 2, geradores diagonais 1 e linear em k). Ele é
repetido aqui como Parte 1, com as variantes abaixo.

## Parte 1: simetria contínua interna da eq. (1)
Álgebra de Lie das transformações ψ → e^{tX}ψ (X anti-hermitiana) que deixam invariante
(ψ̄J^aψ)(ψ̄Σ^aψ) (tensor simetrizado nos dois bilineares), com J^a = J P_L, Σ^a = J P_R
como em WP-J1. Variantes: (a) só rótulos (6×6); (b) rótulos ⊗ quiralidade de Dirac
(X ⊗ 1 e X ⊗ γ5, com o termo G_UV(ψ̄ψ)² também exigido invariante).
Saída: dimensão, geradores, se a álgebra é abeliana.

## Parte 2: quebra pelo condensado de enrolamento (J.6.2)
O condensado gera J⁺Σ⁺ + J⁻Σ⁻. Quais geradores da Parte 1 ele preserva? Variedade de
vácuo G/H e seus π₀…π₃ (para grupos U(1)ⁿ e discretos: π₁ = Zᵏ, π₂ = π₃ = 0).
Critério: existe carga topológica pontual (π₂ ou π₃ ≠ 0)? Se não, (i) está excluída.

## Parte 3: as cargas de Noether como B
Avaliar cada gerador preservado nos estados observados (p, n, Λ, Σ⁰, e, μ, ν's) e
testar se B está no espaço gerado por {geradores, Q, N_ψ}. Se não, (ii) está excluída.

## Parte 4: o nível como grau de liberdade
Contar os estados internos de N = 3 (tuplas de rótulos distintos) e as neutras;
verificar se a eq. (1) tem algum grau de liberdade interno que distinga n de ν_τ
(mesma tupla). Se não, o nível só pode ser excitação espacial, e saltos de nível na
mesma tupla são transições radiativas sem regra de seleção interna: conflito direto
com o setor dual de §6 (n/ν_τ) pela exigência de 5e-30 de WP-LEV1.

## Antecipação registrada
Parte 1: (a) U(1)_N × U(1)_M (M = Σ m_k), abeliana; (b) sem extensão contínua quiral.
Parte 2: o condensado quebra U(1)_M; sobra U(1)_N; π₂ = π₃ = 0: sem skyrmion.
Parte 3: nos estados observados M = 3Q − N_ψ/2 (redundante com Q); B fora do espaço.
Parte 4: só 3 tuplas neutras (as LLR), nenhum rótulo interno de nível: o setor dual é
incompatível com estabilidade, a menos que se acrescente um ingrediente à eq. (1).
