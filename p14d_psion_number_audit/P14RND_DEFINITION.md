# P14 / WP-RND: teste de acaso no dominio continuo do DYN12 (pre-registro)
DECLARACAO: registrado DEPOIS de WP-EXH (que achou ponto viavel em raio, com massas erradas). Objetivo: medir quanto valor tem esse achado.
## Pipeline (identico para alvo real e alvos aleatorios)
(a) Tabela: 12000 pontos (hipercubo latino, mesmo dominio e amostragem de WP-EXH: t,s em [0,6] com 5% em zero e resto log-uniforme em [0,01;6], kappa,lambda em [0,20])
com R_n, R_S0, R_L, R_Sp1, R_Sp2 calculados uma vez (modelo de WP-DYN12/EXH, estado fundamental).
(b) Para um conjunto de alvos (cn,cs0,cl,csp): L = max(|R_n-cn|/0,0164, |R_S0-cs0|/0,15, |R_L-cl|/0,10, min(|R_Sp1-csp|,|R_Sp2-csp|)/0,175).
 Pega-se os 2 melhores pontos da tabela e refina-se cada um por Nelder-Mead (maxiter 60) sobre (t,s,kappa,lambda) recalculando R.
 Sucesso <=> min L < 1.
(c) Alvo real: (-0,164; 0; 0,15; 0,975). Alvos aleatorios: 80 conjuntos com cada centro U[-0,6; 2,2] (mesmas janelas).
## Regra de interpretacao (fixada antes)
P_acaso = fracao de alvos aleatorios com sucesso (IC binomial 95% reportado). Se P_acaso > 0,05 o ajuste de raios do DYN12 NAO e evidencia
(o modelo ajusta alvos arbitrarios com a mesma facilidade). Se <= 0,05, e evidencia fraca de que os raios medidos sao especiais, ainda sem massas.
Limitacao declarada: otimizacao truncada (2 partidas, 60 iteracoes) subestima o sucesso real tanto do alvo real quanto dos aleatorios.
## Antecipacao
P_acaso entre 0,2 e 0,6 (4 parametros continuos contra 4 alvos com janelas largas, exceto n); o alvo real tera sucesso.
