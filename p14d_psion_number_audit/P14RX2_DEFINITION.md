# P14 / WP-RX2: o mesmo mecanismo em 12 posicoes (pre-registro)
Motivo: WP-RX mostrou que com 6 vertices nenhuma tripla Lambda alcanca (0; 0,25] (so R = 0 exato), e 0 de 14641 pontos
(kappa, lambda) cumprem os alvos. Com meia-passos de 30 graus (rotacao fina + troca de psions) WP-GEO9 achou colocacoes
que cumprem todos os raios; a pergunta agora e se UMA funcao de energia as SELECIONA.
## Mecanismo (identico a WP-RX, sitios 12 em vez de 6)
Sitios s = 0..11 no angulo 30 s. Colocacao = 3 sitios distintos para os 3 rotulos (1320). Vertice natural do rotulo k: sitio 2k.
E = E_clock + kappa E_C + lambda E_JS, com
 E_clock = - sum_i cos(30 (s_i - 2 k_i) graus)   (custo de girar/deslocar em meio-passos),
 E_C = + sum_{i<j} Q_i Q_j / d_ij,  E_JS = - sum_{i em L, j em R} 1/d_ij  (d = corda no circulo unitario de 12 sitios).
Grade: kappa, lambda em 0..30, passo 0,5 (3721 pontos). Configuracao do hadron = minimo (empates com raios diferentes -> excluido).
## Observaveis e alvos (iguais a WP-RX)
R_t = <r^2>_t / <r^2>_p, ambos nos seus minimos. n -0,164 a 10%; Sigma0 |R| <= 0,15 (e 0 numerico e permitido: |R| < 1e-9 conta como 0);
Lambda: 0,05 <= R <= 0,25 (fixado agora, para evitar o zero numerico: cobre o intervalo teorico 0,055-0,153 e a margem);
Sigma+ [0,80; 1,15] em (1,3,5) OU (2,3,4). Massas: E_n - E_p > 0; eps fixado por Lambda - p = 177,4 MeV; prever Sigma+ - p (251 MeV)
e Sigma0 - Lambda (77 MeV) a 30%.
Teste de acaso: 2000 conjuntos de alvos aleatorios (centros U[-0,6; 2,2], mesmas meias-larguras: n 0,0164; Sigma0 0,15; Lambda 0,10;
Sigma+ 0,175).
Graus de liberdade: kappa, lambda, escala eps (3); restricoes: 4 raios + sinal + 2 massas (7).
Controle: kappa = lambda = 0 -> colocacao natural (sitios pares), R_n = -0,154.
## Antecipacao (palpite, nao calculado)
Regiao viavel em raios pequena ou vazia (o relogio exige o psion 3 em 180 graus); se existir, as massas nao fecham; P_acaso moderado.
