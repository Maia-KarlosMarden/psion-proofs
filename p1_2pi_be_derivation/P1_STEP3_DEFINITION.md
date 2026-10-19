# P1 — Passo 3: definição operacional de Z_ren (pré-registro)

Registrada em 19/10/2026, antes de qualquer execução do passo 4.
Esta definição não pode ser alterada depois de ver o resultado; mudanças entram
como nova entrada (D3′, D3″ …), com data e motivo, sem apagar esta.

## 1. Ingredientes (todos já definidos nos passos 1 e 2)

- Vértice D1 (geradores completos J^a e Σ^a; J^a P_L e Σ^a P_R; D1-i e D1-ii fixadas).
- Contrações do passo 2: H_S, F_S (termo escalar) e H_J, F_J (termo J·Σ),
  avaliadas com o propagador S em uso.
- Acoplamento comum G para os dois termos, sinal relativo σ ∈ {+1, −1}
  e fator de Dirac N_D ∈ {2, 4}.

## 2. Grade de ramos (pré-registrada)

Os quatro ramos (σ, N_D) ∈ {+1, −1} × {2, 4} são calculados e reportados.
Nenhum é descartado antes do resultado. Cada ramo conta como uma tentativa:
ao comparar com o acaso, o número de tentativas é 4 (ou o produto, se houver
outras escolhas abertas).

## 3. Definição

**(a) Peso da equação de gap.** Sobre o propagador S, o coeficiente do tadpole
na equação de gap 2PI, Σ = δΦ₂/δS, normalizado por N_D², é

    W(S; σ, N_D) = [ N_D² (H_S + σ H_J) − N_D (F_S + σ F_J) ] / N_D²

(Hartree com duas voltas fechadas, +1; Fock com uma volta, −1).

**(b) Propagador crítico.** S_c é a solução auto-consistente da equação de gap no
acoplamento crítico G_c (o menor G com solução não trivial de gap positivo).
W_crit ≡ W(S_c; σ, N_D). Em um primeiro nível (peso unitário e S = projetor
sobre os estados livres) W é computado exatamente; o valor auto-consistente
substitui esse no passo 4.

**(c) Numerador.** N_Fock ≡ 2^(n_livres), com n_livres = 2j o número de estados
livres (k ≥ 1). Isto é uma entrada pré-registrada, não demonstrada; a derivação
é o passo 5. Se o passo 5 falhar, todo resultado que use N_Fock fica marcado
"numerador não provado".

**(d) Fator de renormalização.**

    Z_ren(j; σ, N_D) ≡ N_Fock / W_crit

**(e) Resultado reportado.**

    b_e^(2PI)(j; σ, N_D) ≡ b0 / ( Z_ren · N_c ),   b0 = (4j+3)/(2j+1),  N_c = 2j+1

Só depois de calculado, compara-se b_e^(2PI) com 325/1152 e Z_ren com 32/25.
Nenhum passo do cálculo usa 325/1152, 32/25, m_e ou ρ_Λ.

## 4. Testes de consistência (obrigatórios)

1. Independência de esquema: W_crit com corte (cutoff) e com regularização
   dimensional d = 4 − 2ε; as duas devem coincidir.
2. Família em j: o mesmo procedimento para j = 1/2, 3/2, 5/2, 7/2, 9/2.
3. Sensibilidade à âncora: repetir com k = 0 incluído entre os estados livres.
4. Conservação de carga e simetria Z6 do vértice (verificação numérica).

## 5. Regra de escolha de ramo (pré-registrada)

Um ramo só pode ser adotado por critérios independentes de ele reproduzir
32/25: (i) existência de solução de gap positiva e G_c > 0; (ii) independência de
esquema; (iii) conservação de carga e simetria Z6; (iv) coerência com outros
resultados registrados da PH. Se mais de um ramo passar (i) a (iv), todos são
reportados e nenhum é adotado. Escolher o ramo que dá 32/25 é circular e conta
como falha de auditoria.

## 6. Critérios de resultado

Valem os da seção 5 do documento de especificação. Em particular, um W_crit
diferente de (2j)² = 25 sem forma fechada em j é registrado como achado negativo.

## 7. Aviso de transparência

Ao redigir esta definição foram calculados, antes do registro, os valores no
limite de peso unitário (S = projetor sobre os livres), apenas como aritmética
dos números do passo 2:

| σ | N_D | W (peso unitário) |
| --- | --- | --- |
| +1 | 2 | 45/4 |
| +1 | 4 | 109/8 |
| −1 | 2 | 135/4 |
| −1 | 4 | 271/8 |

Nenhum é igual a 25. Esses valores não são o resultado do passo 4, porque o
propagador auto-consistente tem pesos diferentes de 1; mas indicam que, para
a definição acima reproduzir (2j)², é preciso que o propagador crítico ou a
dominância do figura-8 mudem os pesos de forma não trivial.
