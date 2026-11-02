# P14 / WP-SEL: criterio de selecao por media termica (Boltzmann) sobre as 120 colocacoes (pre-registro)
DECLARACAO: registrado DEPOIS de RX, RX2, RX3 e EXC (todos falharam com "minimo de energia"). Nao e cego quanto a Lambda ser o bloqueio.
## Criterio
Hexagono de 6 vertices (vertice s a 60 s graus). Cada tripla ocupa 3 dos 6 vertices: 120 colocacoes (permutacoes). Energia como em WP-RX:
E = -sum cos(theta_vertice - theta_natural) + kappa sum q_a q_b/d_ab - lambda sum_{pares L-R} 1/d_ab.
Em vez do minimo: <R2> = sum_i w_i R2_i / sum_i w_i, w_i = exp(-(E_i - E_min)/T). R = <R2>_h / <R2>_p (mesma T, kappa, lambda).
Grade: kappa, lambda em 0..30 (passo 1,0); T em {0,1; 0,3; 1; 3; 10; infinito(media uniforme)}: 31*31*6 = 5766 pontos. Parametros livres: kappa, lambda, T (3).
Alvos e janelas: n -0,164 +-0,0164; Sigma0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] em (1,3,5) ou (2,3,4).
Teste de acaso: 2000 conjuntos de alvos aleatorios (centros U[-0,6;2,2]).
## Antecipacao
Raios: 0 pontos viaveis (Lambda ou Sigma0 deixa de ser 0 fora do minimo; media uniforme nao respeita os alvos); P_acaso baixo mas > 0 por causa de T.
Se houver pontos viaveis, serao em T intermediario, e seriam o primeiro avanco desde WP-HUB; massas nao testadas neste WP.
