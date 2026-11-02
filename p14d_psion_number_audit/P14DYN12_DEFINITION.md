# P14 / WP-DYN12: dodecagono dinamico (estado fundamental quantico) (pre-registro)
DECLARACAO: registrado DEPOIS de RX, RX2, RX3, EXC, SEL (todos falharam em Lambda). Estende WP-DYN (hexagono, 120 estados) a 12 sitios.
O dodecagono foi originalmente proposto para neutrinos; aqui e aplicado a barions como teste da geometria mais fina. Nao e cego.
## Modelo
Estados: colocacoes ordenadas dos 3 rotulos da tripla em 3 dos 12 sitios (1320). Sitio s a 30 s graus; sitio natural do rotulo k = 2k.
H = diag(E_clock + kappa E_C + lambda E_JS) + hopping:
 -t : um psion salta +-1 sitio (30 graus) ou +-2 sitios (60 graus, angulo hexagonal) para sitio livre;
 -s : troca de posicao entre dois psions OU salto de 6 sitios (180 graus) para sitio livre (mecanismo antipodal).
E_clock = -eps sum cos(theta_sitio - theta_natural) com eps=1; E_C, E_JS como em WP-RX (B = sum q q/d; C = -sum_{L-R} 1/d).
Estado fundamental (menor autovalor, via eigsh); R2_h = sum |psi_i|^2 R2_i; R = R2_h/R2_p (mesmos parametros).
Grade: t, s em {0; 0,1; 0,3; 1; 3}; kappa, lambda em {0; 2; 8}: 225 pontos. Parametros livres: t, s, kappa, lambda.
Alvos/janelas: n -0,164 +-0,0164; Sigma0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] em (1,3,5) ou (2,3,4). Teste de acaso: 2000 alvos aleatorios.
Degenerescencia exata do fundamental (t=s=0, empates): resultado marcado NaN (excluido), como em RX.
## Antecipacao
0 pontos viaveis; Lambda continua bloqueando (com t pequeno, localizado, R_L fica em valores discretos do arco; com t grande, tende ao valor
deslocalizado, onde n falha como em WP-DYN). Se houver ponto viavel, sera em regime intermediario e primeiro avanco desde WP-HUB.
