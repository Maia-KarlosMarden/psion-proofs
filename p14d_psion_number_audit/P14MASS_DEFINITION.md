# P14 / WP-MASS: termos de massa com coeficientes livres + massa antes dos raios (DYN12) (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: motivado por WP-MFIRST (kappa>=0 da m_n<m_p e m_L<m_p em todos os pontos). Registrado depois dele. Nao cego.
## Termos acrescentados (coeficientes livres)
E_h = E_gs(t,s,kappa,lambda) + mu * M_h + nu * Q_h, com M_h = sum (k - 5/2) sobre os rotulos (p: +1,5; n, S0, L: -1,5; Sp1, Sp2: +1,5) e
Q_h = carga total (p: 1; n, S0, L: 0; Sp: 1). Sao constantes por hadron: NAO alteram funcao de onda nem raios.
Consequencia exata: so a combinacao c = -(3 mu + nu) importa: desloca {n,S0,L} em c relativo a {p,Sp}. Portanto:
 E_n-E_p = dgs_np + c ; E_L-E_p = dgs_Lp + c ; E_Sp-E_p = dgs_Spp ; E_S0-E_L = dgs_S0L ; E_L-E_n = dgs_Ln (todas inalteradas por c exceto as duas primeiras).
## Criterios (massa PRIMEIRO)
D = E_gs_L - E_gs_n > 0 (so o estado fundamental; c nao ajuda). Escala u = 176,1 MeV / D (Lambda - n = 1115,68 - 939,57 = 176,1).
R1: (E_S0-E_L)/D = 77/176,1 = 0,437 +-30%. R2: (E_Sp-E_p)/D = 251/176,1 = 1,425 +-30% (Sp1 ou Sp2).
c e entao FIXADO por n-p = +1,29 MeV: c = 1,29/u - (E_gs_n - E_gs_p); reporta-se |c| u em MeV e |c|/D (naturalidade: |c|/D <= 3 considerado nao ajustado fino).
Depois: entre os pontos que passam D>0, R1, R2, o menor L de raios (perda de WP-EXH; viavel <=> L<1) e a contagem de L<1 e L<2.
Dominio fisico: kappa em [0,20], lambda em [0,20], t,s em [0,6] (5% em zero, resto log-uniforme [0,01;6]), 12000 pontos LHS.
Controle de rigidez: fracao de pontos que passam R1&R2&D>0 (e fracao de alvos aleatorios de massa nao testada).
## Antecipacao
R1&R2 e D>0 passam em poucos pontos (<2%); entre eles L de raios > 1 (0 pontos com L<1). Se houver, c sera grande (|c|/D > 3).
