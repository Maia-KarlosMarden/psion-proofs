# P14 / WP-C: sphalerons na temperatura de nucleação (pré-registro)

Registrado antes de rodar. Entradas: Λmeta = 1.2609e12 GeV, ξ* = 434.09 (entrada
PROVISÓRIA, só para T_nuc = Λmeta/ξ*; será substituída por WP-A), g* = 106.75,
M_Pl = 1.2209e19 GeV, α_2(M_Z) = 0.0338 com corrida a 1 laço do SM (b_2 = 19/6),
taxa de sphaleron Γ/T⁴ = κ α_W⁵ com κ = 18 ± 3 (fase simétrica, valor de
referência da literatura de lattice; registrar a fonte na redação final).

## Passo 1 (leitura, feito antes deste registro)
O paper não declara T_gauge. §18.2 afirma a derivação do grupo de gauge sem
escala; Definição I.1 confina proto-quarks abaixo de Λmeta; §13.1 diz que
SU(2)meta "não tem relação com o SU(2) eletrofraco". Logo os casos (a)/(b)/(c)
não são decidíveis pelo texto. Resultado do passo 1: INDETERMINADO; a condição
T_gauge > T_nuc entra como premissa declarada (P14-b).

## Passo 2
Calcular Γ_sph/(T³H) em T_nuc e em Λmeta, com α_W corrido. Critério: > 1 ativo.

## Passo 3
Tentar reproduzir Ccore = 8.4434e-4 a partir da fórmula impressa em §13.3,
η_sph = (8/23)(g*/g_S^{4/3})×1e-4, com g_S ∈ {g*(4/11), 3.91, g*, 43/11·...}
e variações óbvias. Se nenhuma reproduzir a 1%, registrar "não reprodutível".
Verificar separadamente se a eq. (14) com fatores O(1) pode dar 6.02e-10.

## Passo 4
Diluição de entropia: n_B/n_γ em T_nuc vs hoje, fator g_S(T_nuc)/g_S(hoje).
Verificar se a eq. (33) contém esse fator (n/n_γ = 0.749 ≈ 3/4 sugere razão
férmion/bóson por grau de liberdade, sem diluição).

## Passo 5
η_B pela cadeia de WP-B (c_s(B−L) no lugar de ε_conf) com ξ* = 434 provisório,
para os esquemas S1 e S3, com e sem o fator 1/3 por psion, com e sem diluição.
Só comparação final com 6.13e-10.

## Critério pré-registrado
Sucesso parcial é o máximo possível (passo 1 indeterminado). Reportar Ccore
reprodutível ou não; ativo ou não em T_nuc e em Λmeta.
