# Parecer de referee: bariogênese na PH (WP-BIF a WP-COMP, artigo e §13 do v86)

Autor do parecer: o mesmo assistente que construiu os testes (não é revisão independente).
Data: 20/10/2026. Checagem numérica exploratória deste parecer: `P14REF_check_dl.txt`.

## Recomendação
**Revisão maior.** Os resultados estruturais são sólidos e publicáveis. A parte dinâmica,
e em especial a "previsão" T_RH* ≈ (1–2)×10¹⁰ GeV, não se sustenta como está escrita.

## O que está correto e é forte
1. Resultados exatos: conservação de N_ψ, simetria U(1)_N × U(1)_M, CP de campo exata,
   ontologia O−, B e L exigidos com N_ψ/3 = B + L, ausência de B de Noether ou topológico.
2. Estrutura da bifurcação: B = −L, Q = 0 e o estado final hidrogênio + 2ν̄ decorrem de uma
   simetria das equações, com prova e verificação numérica.
3. Seleção: as 8 triplas singleto de proto-cor coincidem exatamente com o critério de
   confinamento do próprio paper (I mod 6 ∈ {0,3}). As exóticas não são problema.
4. Disciplina: pré-registro com hash, falhas e erros documentados.

## Problemas maiores
**M1. Atribuição de degraus inconsistente entre testes.** WP-UNW exigiu que as cópias
leptônicas carregadas positivas sejam mais pesadas que o próton (ℓ_L ≤ 7). WP-RUNG e
WP-COMP usaram ℓ_L = 10 para todas as carregadas (Δℓ = +2). Com a atribuição consistente
(Δℓ = +2 só para a tupla do elétron, +7 para as neutras e −1 ou −2 para as positivas), os
pesos das positivas trocam de sinal. A checagem deste parecer mostra que o sinal de η passa
a mudar com T_RH e com a lei de ligação: em L1, −7,5e-5 em 7e8 e +7,7e-3 em 7e9; em L2,
+9,8e-5 em 7e9 e −2,8e-6 em 1,4e10. Como o sinal global já é um acaso do vácuo (k ou 6−k),
isso não mata o mecanismo. Mas |η| passa por zero em várias temperaturas, então surgem
raízes T_RH* adicionais, e a afirmação de uma faixa única e robusta cai.

**M2. Dependência da condição inicial.** Em WP-COMP a assimetria nasce só porque a rede parte
de zero compósitos em T_RH (reaquecimento instantâneo em radiação pura). Se o ínflaton decai
gradualmente (com temperatura máxima T_max > T_RH), ou decai em parte direto em compósitos,
a fonte e a borda de lavagem mudam: a borda passa a ser controlada por T_max, não por T_RH.
O modelo de reaquecimento é parte do mecanismo e não foi modelado.

**M3. Normalização das taxas.** T_RH* fica na borda de lavagem, fixada por T_eq ∝ c_Γ^{−1/3}.
O valor c_Γ = 1 ("acoplamento natural") supõe que a criação de pares de compósitos segue a
taxa de contato G4 T⁵/Λ⁴ com G4 = (2π²/3)². Um fator 10 em c_Γ desloca T_RH* por ~2,2. A
"robustez" vale para ε, não para a normalização da taxa.

**M4. Sintonia escondida.** No pico, o mecanismo superproduz: |η| ~ 1e-2 contra 6,1e-10,
um excesso de 10⁷. Ele só reproduz o valor observado se (a) ε ~ 1e-8 (sem explicação) ou
(b) T_RH estiver a ~±15% da borda de lavagem. A opção (b) é a "previsão" de T_RH, mas também
é uma sintonia: o universo precisaria ter reaquecido exatamente na borda. Um referee externo
diria isso de imediato.

**M5. Fonte de CP com duas ambiguidades.** (i) Em uma das três colocações da fase do vácuo
(PA) ela é removível e o efeito some; o paper não fixa a colocação. (ii) O vácuo de menor
energia (s = 0) é invariante sob CP, e o mecanismo exige que a região inflacionária esteja
num vácuo metaestável k ∈ {1, 2, 4, 5}, um acaso cosmológico. Além disso, a assimetria
dinâmica de WP-SGN foi medida na rede com psions livres, que M7 invalida; só a parte
algébrica (fase física, ímpar, nula em k = 0, 3) sobrevive.

**M6. Caminho adaptativo.** Cada etapa foi pré-registrada, mas a sequência foi adaptativa:
calibração → bifurcação → T_RH → sinal → degrau → compósitos, cada hipótese formulada depois
de ver o resultado anterior. O pré-registro por etapa não protege contra a escolha sequencial
do sobrevivente. Falta uma previsão fora da amostra, isto é, algo calculado antes de comparar
com qualquer dado.

**M7. Resultados herdados da rede com psions livres.** WP-SIM, NOW, DOM, CAL, OPS, LOOP, BIF e
TRH usaram psions livres abaixo de Λ_meta, o que o App. I proíbe. As conclusões estruturais
sobrevivem (são simetrias), mas os números não: a tabela de eficiência do artigo e do §13.7
do v86 vem de WP-TRH e precisa ser substituída pela de WP-COMP, com as ressalvas M1–M4.

## Problemas menores
- g* = 106,75 (Modelo Padrão) usado a 10¹⁰ GeV numa teoria com outro espectro; efeito
  ~g*^{1/6} em T_eq, pequeno, mas deve ser declarado.
- Cada degrau tratado como espécie com g = 2: dobra os graus de liberdade das tuplas.
- Limiares absolutos de controle (1e-12, 1e-18) geraram falhas formais no piso numérico.
  Registrar limiares relativos ao sinal.
- κ (parte absortiva) escolhido entre três valores "naturais" sem derivação.
- O artigo cita T_RH ≈ 10⁹–4×10¹⁰ GeV como previsão (5); com M1–M4 isso deve virar
  "condição de consistência", não previsão.

## Testes pedidos antes de publicar a parte dinâmica
1. Refazer WP-COMP com a atribuição consistente de M1 e listar todas as raízes T_RH*.
2. Modelar o reaquecimento (decaimento do ínflaton com Γ_φ, T_max > T_RH) e medir a
   dependência de T_max.
3. Varrer c_Γ em uma década e reportar T_RH*(c_Γ).
4. Fixar a colocação da fase do vácuo a partir do próprio J.6.2, ou declarar a ambiguidade
   como parâmetro.
5. Formular uma previsão fora da amostra antes de qualquer nova comparação.

## Avaliação
A parte estrutural e os no-go são um artigo publicável hoje. A parte dinâmica é um mecanismo
candidato coerente com a PH, mas com sintonia escondida (M4) e sem previsão firme (M1–M3).
Placar revisado da bariogênese: ~55%, não 70%.
