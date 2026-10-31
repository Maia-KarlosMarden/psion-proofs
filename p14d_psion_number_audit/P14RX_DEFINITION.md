# P14 / WP-RX: rotacao do hexagono + troca de psions como mecanismo de selecao de configuracao (pre-registro)
Pedido do autor: conceber um mecanismo de rotacao do hexagono combinado com troca de psions que ajuste, de uma so vez,
os valores que destoam (raios de carga e massas).
## Mecanismo
Cada tripla t ocupa uma das 120 colocacoes dos seus 3 rotulos nos 6 vertices (rotacao r do hexagono + permutacao/troca
de posicoes, incluindo reflexao). O hadron e a colocacao de MENOR energia (configuracao classica; sem tunelamento).
E(colocacao) = E_clock + kappa E_C + lambda E_JS, em unidades de eps = 1:
 E_clock = - sum_i cos(theta_{v_i} - theta_{k_i})      (relogio Z6: custo de girar/trocar em relacao ao vertice natural;
                                                         rotacao rigida de r passos custa 3(1 - cos 60 r))
 E_C     = + sum_{i<j} Q_i Q_j / d_ij                  (eletrostatica de cordas do anel, raio 1)
 E_JS    = - sum_{i em L, j em R} 1 / d_ij             (acoplamento L-R do termo J.Sigma, atrativo)
Parametros livres: kappa >= 0 e lambda >= 0 (grade 0..30, passo 0,25). A escala eps em MeV e fixada por (E_Lambda - E_p).
## Observaveis
Raios: R_t = <r^2>_t(min) / <r^2>_p(min) (o proton tambem escolhe sua colocacao). Alvos de WP-GEO8/9: n -0,164 a 10%;
 Sigma0 |R| <= 0,15; Lambda (0; 0,25]; Sigma+ [0,80; 1,15] em (1,3,5) OU (2,3,4); Sigma- 0,86 +- 0,21 (opcional: identificacao).
Massas (so diferencas): E_n - E_p > 0 (sinal); apos fixar eps por Lambda - p = 177,4 MeV, prever
 Sigma+ - p (251 MeV) e Sigma0 - Lambda (77 MeV), a 30%. (Sigma- fora: mesma tripla do eletron.)
## Contagem de graus de liberdade
2 parametros de forma (kappa, lambda) + 1 escala (eps) = 3 livres. Restricoes: 4 raios (n, Sigma0, Lambda, Sigma+) + 1 sinal + 2 massas = 7.
## Teste de acaso (look-elsewhere), registrado antes
Repetir a busca com alvos ALEATORIOS: para cada um dos 4 estados, centro sorteado em U[-0,6; 2,2] com as MESMAS meias-larguras
(n 0,0164; Sigma0 0,15; Lambda 0,125; Sigma+ 0,175). 2000 sorteios. P_acaso = fracao com algum ponto da grade viavel.
## Controle
kappa = lambda = 0: colocacao natural; R_n = -0,154, R_Sigma- = 0,846 (hexagono estatico).
## Antecipacao (palpite, nao calculado)
Existe regiao viavel nos raios (varios pontos); P_acaso e alta (>30%); as massas nao fecham (Sigma+ - p e Sigma0 - Lambda fora de 30%).
