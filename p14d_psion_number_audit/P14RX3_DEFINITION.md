# P14 / WP-RX3: o mecanismo de 12 posicoes com o vinculo L-R 2-3 (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: este termo foi escolhido DEPOIS de ver que o Lambda e o bloqueio de WP-RX/RX2 e de saber
que a solucao de WP-GEO9 poe o psion 3 a 150 graus, junto do 2. Portanto e um ajuste dirigido, nao um teste cego.
Motivacao independente (so para registro): o App. J.6.2 do paper gera um vinculo real |2>-|3> (h1[2,3] = 0,5).
## Mecanismo
Igual a WP-RX2 (12 sitios, vertice natural do rotulo k = sitio 2k, E = E_clock + kappa E_C + lambda E_JS) mais
 E_bond = - beta / d_23   (atracao entre os psions de rotulos 2 e 3, so em triplas que contem ambos: Lambda (1,2,3)
 e Sigma+ (2,3,4)).
Grade: kappa, lambda em 0..30 (passo 0,5), beta em 0..30 (passo 0,5): 226.981 pontos.
## Alvos e criterios: os de WP-RX2 (n -0,164 a 10%; Sigma0 |R| <= 0,15; Lambda [0,05; 0,25]; Sigma+ [0,80; 1,15] em (1,3,5) ou (2,3,4)).
Massas: E_n - E_p > 0; eps fixado por Lambda - p = 177,4 MeV; prever Sigma+ - p (251) e Sigma0 - Lambda (77) a 30%.
Graus de liberdade: kappa, lambda, beta, eps (4); restricoes 7.
Teste de acaso: 2000 conjuntos de alvos aleatorios (centros U[-0,6; 2,2]; meias-larguras 0,0164; 0,15; 0,10; 0,175).
## Antecipacao (por argumento, antes de rodar)
Raios: regiao viavel existe em beta moderado (a atracao puxa o 3 para o 2 so no Lambda e no Sigma+(2,3,4), sem tocar p, n, Sigma0).
Massas: FALHAM: o vinculo (e a eletrostatica com kappa > 0, que liga mais as triplas neutras e negativas) abaixam E_Lambda e E_Sigma+
abaixo de E_p; com kappa = lambda = 0 todas as energias naturais sao iguais (E_n - E_p = 0, sem sinal). P_acaso sobe em relacao a RX2.
