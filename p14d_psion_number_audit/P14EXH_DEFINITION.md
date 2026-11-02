# P14 / WP-EXH: esgotamento do dodecagono dinamico (DYN12) (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: motivada por WP-SWP (n+Lambda juntos em 4 pontos; Sigma+ nunca). Busca dirigida a achar QUALQUER ponto viavel.
## Modelo
Identico a WP-DYN12/SWP (12 sitios, 1320 colocacoes, saltos +-1/+-2 com -t, troca/salto antipodal com -s, eps=1, estado fundamental).
## Funcao-perda
d_n=|R_n+0,164|/0,0164; d_S0=|R_S0|/0,15; d_L=|R_L-0,15|/0,10; d_Sp=min(|R_Sp1-0,975|,|R_Sp2-0,975|)/0,175; L=max(d_n,d_S0,d_L,d_Sp).
Viavel <=> L<1. NaN (degenerescencia) -> L=50.
## Dominio (continuo)
t, s em [0, 6] (amostragem: 0 com prob. 5%, senao log-uniforme em [0,01; 6]); kappa, lambda em [0, 20] (uniforme).
## Procedimento
(1) Hipercubo latino com 6000 pontos. (2) Evolucao diferencial (scipy), 3 sementes, popsize 12, maxiter 60, sobre (t, s, kappa, lambda) em [0,6]x[0,6]x[0,20]x[0,20].
(3) Nelder-Mead (maxiter 200) a partir dos 8 melhores pontos de (1) e dos 4 pontos n+Lambda de WP-SWP (reencontrados: pontos da grade SWP com n e Lambda ok; se nao reproduzidos, dos 8 melhores). Imprimir R_n, R_S0, R_L, R_Sp1, R_Sp2 de cada otimo.
## Criterio de "esgotado" (fixado antes)
Se min L >= 1 em todas as etapas, declara-se: NENHUM ponto viavel encontrado em [0,6]^2x[0,20]^2 por busca global + local; o DYN12 esta esgotado como mecanismo
dentro deste dominio (nao prova de ausencia global). Se L<1 for achado: checar massas (criterios de WP-SWP) e P_acaso(=0,003 da grade densa nao se aplica; novo teste de acaso nao e feito).
## Antecipacao
min L > 1 (provavelmente entre 1,2 e 3, limitado por Sigma+ ou n); nenhum ponto viavel.
