# P1 — D3′: redefinição operacional de Z_ren (pré-registro)

Entrada nova no registro de decisões. Não apaga nem altera D3
(P1_STEP3_DEFINITION.md, SHA-256 4e0ad99b…2a2c). D3 e seu resultado (passo 4)
permanecem válidos como registro; D3′ é a segunda tentativa.

Motivo: o passo 4 expôs um defeito de D3, a normalização D3-N (Tr s_c = n_livres).
W é forma quadrática em s, então W[s_c] depende da escala e da forma da
normalização; no ramo σ = +1 o autovetor crítico tem traço nulo e a
normalização fica indefinida. A escala tem de ser fixada pela equação de gap,
não por convenção.

## 1. Análise (antes de definir)

1. W[s] = ⟨s, 𝒲 s⟩, com ⟨a,b⟩ = Tr(a b) (produto de Frobenius, o mesmo
   emparelhamento que aparece na ação efetiva Tr(M s)). 𝒲 é um operador linear
   simétrico fixo no espaço das matrizes Hermitianas de rótulos.
2. Equação de gap não linear: M = g 𝒲[s], s = M r(M²), com r(x) = f(x)/f(0) e
   f a integral de laço regularizada (corte: f ∝ Λ² − x ln(1+Λ²/x); tempo
   próprio: f ∝ Λ² e^(−x/Λ²) − x Γ(0, x/Λ²)); g = 2 G f(0) / (normalização
   de Dirac).
3. Linearizando: g_c = 1/λ_max(𝒲). O λ_max é o autovalor do operador 𝒲 no
   emparelhamento de Frobenius. Ele **não depende da normalização** do
   autovetor, e não depende do regulador (que só entra em r e em g_c).
4. Logo acima de g_c a solução não linear é s* = μ v + O(μ³), com v o autovetor
   crítico (unitário em Frobenius) e μ² ln(1/μ²) ∝ (g − g_c)/g_c (corte).
   Φ₂(g) ∝ μ² → 0, e d²Φ₂/d(ln g)² também → 0 quando g → g_c⁺. Portanto
   a segunda derivada pura de Φ₂ não define um número finito. O que tem limite
   finito e independe de escala é a razão W[s*]/Tr(s*²) → λ_max.

## 2. Definição D3′

    W_crit′(j; σ, N_D) ≡ n_livres · lim_{g→g_c⁺} W[s*(g)] / Tr(s*²)
                       = n_livres · λ_max(𝒲; σ, N_D)

    Z_ren′ ≡ 2^(n_livres) / W_crit′
    b_e′   ≡ b0 / (Z_ren′ · N_c),  b0 = (4j+3)/(2j+1),  N_c = 2j+1

O fator n_livres é uma escolha de unidade, declarada: Tr s² = n_livres é a
norma de Frobenius do projetor de peso unitário sobre os estados livres, que é
a referência do passo 2 (para ele W = n² − n/N_D no termo escalar). Com isso
D3′ coincide com D3 no limite em que a direção crítica é o projetor unitário.
Não será tentada outra unidade sem nova entrada (D3″).

D1, D1-i … D1-iv, D4-G e a grade de quatro ramos (σ, N_D) ficam como estavam.

## 3. Grade e contagem de tentativas

Os quatro ramos são reportados, nenhum adotado (regra da seção 5 de D3).
Tentativas acumuladas para correção de busca: D3 (4) + D3′ (4) = 8, mais as
variantes já reportadas (4 com carga restrita, 2 com âncora). Para efeito de
look-elsewhere o N usado é 8 nos ramos principais.

## 4. Verificações obrigatórias do passo 4b

1. Resolver a equação de gap NÃO LINEAR (Newton, jacobiano numérico) para
   g = g_c(1+δ), δ ∈ {1e-1, 1e-2, 1e-3, 1e-4}, com dois reguladores (corte e
   tempo próprio), e verificar que W[s*]/Tr(s*²) → λ_max e que s* tem
   sobreposição → 1 com o autoespaço crítico.
2. Verificar numericamente que d²Φ₂/d(ln g)² → 0 (justifica não usá-la).
3. Independência de esquema: a razão-limite coincide nos dois reguladores.
4. Família em j = 1/2 … 9/2.
5. Lado a lado com D3 (W_crit de D3-N).

## 5. Aviso de transparência (resultado já parcialmente conhecido)

Antes de registrar D3′, a tabela de n_livres·λ_max foi calculada a partir das
autovalores já impressos na saída do passo 4 (sem rodar nada novo). Em j = 5/2:

| σ | N_D | λ_max | Z_ren′ = 2^5/(5 λ_max) |
| --- | --- | --- | --- |
| +1 | 2 | 2.9946 | 2.137 |
| +1 | 4 | 3.5510 | 1.802 |
| −1 | 2 | 6.8099 | 0.9398 |
| −1 | 4 | 6.8729 | 0.9312 |

Nenhum é 32/25 = 1.28. Logo D3′ é registrada sabendo que não reproduz Z_NP;
ela é rodada para (a) validar a solução não linear, (b) testar esquema,
(c) obter a família em j. O resultado não é cego e não conta como evidência
a favor, nem contra, de qualquer ramo além do que o número mostra.
