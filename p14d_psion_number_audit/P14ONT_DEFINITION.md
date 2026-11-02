# P14 / WP-ONT: qual antimatéria é consistente com reações observadas (pré-registro)

Pergunta do autor: a antipartícula tem +3 ou −3 psions? Teste: contabilidade exata
sobre um catálogo de reações medidas. A eq. (1) conserva o número de psions N_ψ
exatamente, e a carga elétrica é conservada nas reações observadas.

## Ontologias
- O+ (A.5 lida literalmente): antipartícula = tupla com índices invertidos
  k → (6−k) mod 6, feita de psions: N_ψ = +3 para matéria e antimatéria de três
  corpos; carga = soma de Q9 dos índices invertidos.
- O− (antipsions): antipartícula = mesma tupla feita de antipsions: N_ψ = −3,
  carga = −Q, índice Z6 invertido (ψ̄ carrega Z6 conjugado). Nesta leitura a
  inversão de A.5 descreve o índice Z6 da antipartícula, não o número de psions.

## Atribuições (do paper)
e⁻ = (0,1,2) nível 10; μ⁻ = excitação de (0,1,2); p = (0,4,5) nível 8;
n = (0,1,5) nível 8; ν_e = (1,2,3), ν_μ = (0,2,4), ν_τ = (0,1,5) (§6, setor dual).
Q9 = (1,−5,−5,10,4,4)/9. Bósons (γ, π⁰, π⁺, π⁻, X17) com N_ψ DESCONHECIDO (livre),
carga observada. Modo "fiel ao paper": X17 tem N_ψ = 0 (K.8.2: setor misto ψψ̄).

## Catálogo de reações observadas
R1 n → p e⁻ ν̄_e; R2 p → n e⁺ ν_e (β⁺ em núcleos); R3 e⁻ p → n ν_e (captura);
R4 e⁺e⁻ → γγ (para-positrônio); R5 e⁺e⁻ → γγγ (orto-positrônio);
R6 γ → e⁺e⁻ (conversão perto de um núcleo; núcleo inalterado);
R7 p p̄ → π⁺π⁻π⁰; R8 π⁰ → γγ; R9 π⁺ → μ⁺ν_μ; R10 μ⁻ → e⁻ ν̄_e ν_μ;
R11 X17 → e⁺e⁻ (previsão da própria PH, K.8.2; Atomki).

## Testes
T-N: existe atribuição de N_ψ aos bósons que conserva N_ψ em todas as reações?
     (sistema linear exato; reportar quais reações conflitam)
T-Q: a carga se conserva em todas as reações com as cargas das tuplas?
T-A: propriedades de antipartícula observadas: carga oposta para e, p, Σ⁺;
     n̄ distinto de n (n e n̄ se aniquilam; momentos magnéticos opostos).

## Critério
Uma ontologia é consistente se passa T-N, T-Q e T-A. Reportar as duas. Se só uma
passar, ela é a recomendada para a contabilidade dos testes da regra de níveis.

## Transparência
Ao redigir, antecipei que O+ falha em R4 + R5 (6 psions não se dividem em 2 e em 3
fótons com o mesmo N_ψ por fóton) e que, em O+, o nêutron (0,1,5) seria a própria
antipartícula. O script verifica isso e o resto do catálogo.
