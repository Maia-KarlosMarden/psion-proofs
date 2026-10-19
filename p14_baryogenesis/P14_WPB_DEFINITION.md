# P14 / WP-B: contabilidade de N_ψ, B e L (pré-registro)

Registrado antes de rodar. Entradas proibidas: η_B, ξ*, Ω_DM.

## Achados de leitura (antes de qualquer conta)

1. A Tabela 2 do paper lista 8 compostos observáveis de três corpos
   (I = Σk mod 6 ∈ {0,3} entre as C(6,3) = 20 triplas). A nota de ε_conf fala
   em "nove famílias".
2. ε_conf = 4/9 tem três justificativas distintas no material:
   (a) N3/(N2+N3+N4), com N2, N3, N4 não definidos no texto;
   (b) "quatro compostos leptônicos (e, μ, τ, ν_e) entre nove famílias";
   (c) (2/N)² com N = 3 (P_yield_ratio.py, "confinement scaling").
   Nenhuma é derivada. Isto é registrado como achado, independentemente do
   resultado abaixo.

## Definições

- Censo: regenerar as 20 triplas, filtrar I ∈ {0,3}, conferir que dá a Tabela 2.
- Cargas por psion Q9 = (1,−5,−5,10,4,4)/9 (Teorema V.1, convenção do paper).
- N_ψ(composto) = 3. Assimetria primordial Δn_ψ (psions menos antipsions),
  de sinal fixado por ε_ψ; aqui tratado como unidade: tudo em Δn_ψ/3.
- Classificação (fixada antes de contar). Três esquemas, todos reportados:
  S1 "leitura do paper": 9 famílias = {e, μ, τ, ν_e(1,2,3)} léptons;
     {p(0,4,5), n(0,1,5), Σ⁰(0,2,4), Σ⁺(1,3,5), Σ⁺(2,3,4)} bárions; exótico
     (3,4,5) excluído. Gerações contam como famílias.
  S2 "8 triplas, duais = neutrinos": léptons {(0,1,2),(0,1,5),(0,2,4),(1,2,3)};
     bárions {(0,4,5),(1,3,5),(2,3,4)}; exótico B = L = 0.
  S3 "8 triplas, duais = bárions": léptons {(0,1,2),(1,2,3)}; bárions os outros
     cinco; exótico B = L = 0.
  Peso: democrático f_i = 1/N_fam; e, como sensibilidade, peso 2J+1 com
  J = |Σm_J|.
- B = Σ_bárions f_i, L = Σ_léptons f_i (em unidades de Δn_ψ/3).
- Conversão por sphaleron: B_final = c_s (B − L), c_s = (8N_f + 4N_H)/(22N_f + 13N_H);
  Modelo Padrão N_f = 3, N_H = 1: 28/79. Testar se o 8/23 de §13.3 é c_s para
  algum (N_f, N_H).

## Teste de conservação (passo 2)

Entre todos os pares de compostos observáveis, listar as recombinações 2→2 que
conservam o multiconjunto de k (logo N_ψ e Z6) e verificar se B, L e B−L são
conservados em cada esquema. Hipótese: B e L não são conservados pelo
rearranjo; B−L também não (só N_ψ é).

## Critério (pré-registrado)

Sucesso: passo 2 confirma que o vértice/rearranjo não conserva L, e B−L ≠ 0 com
pesos democráticos no esquema S1. Reportar o fator multiplicativo
c_s (B−L) em relação ao ε_conf = 4/9 da eq. (33), sem ajustar nada. Se o fator
mudar η_B em mais de um fator 3, isso é registrado como mudança da previsão.
