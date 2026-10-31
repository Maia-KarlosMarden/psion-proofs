# P14 / WP-OCT: octogono dinamico (pre-registro)
DECLARACAO: registrado depois do fechamento do DYN12 (P14DYN12_CLOSURE.md); conhece-se a fronteira massa-raio do dodecagono. Pedido do usuario. Nao cego quanto a essa fronteira.
## Modelo
8 sitios (45 graus), antipoda = +4. Rotulos k=0..5 em sitios HOME = (0,1,2,4,5,6): a antipoda de k (k+3) cai em +4; sitios 3 e 7 vagos. L = k<=2 (sitios 0,1,2), R = k>=3 (sitios 4,5,6).
336 colocacoes ordenadas por tripla (8*7*6). H = diag(relogio eps=1: -sum cos(theta_sitio - theta_home(k)); kappa*Coulomb sum q q/d; lambda*J.Sigma = -sum_{pares L-R} 1/d) + saltos:
-t: um psion salta +-1 sitio (45 graus) ou +-2 sitios (90 graus) para sitio livre; -s: troca de posicao entre dois psions OU salto de +4 sitios (180 graus) para sitio livre. Estado fundamental (eigsh).
Raio: <r^2> de carga em torno do centroide / Q_tot (Q_tot=0: sem normalizar), razoes a <r^2>_p. Massa M_h = a<A>+b<B>+g<C> com direcao natural (mesmo sinal).
## Criterios (inalterados de WP-ESC)
Raios: L_r<1 (perda de EXH: n -0,164 +-0,0164; S0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] em (1,3,5) ou (2,3,4)).
Massa: D=M_L-M_n>0; (S0-L)/D em 0,437 +-30%; (Sp-p)/D em 1,425 +-30%; escala natural q=D/max|d_i|>=0,336 (eps<=349 MeV, |w|<=1,5); sensibilidade q>=0,168; 4000 direcoes Dirichlet +/-.
## Procedimento
1. Tabela LHS de 12000 pontos (t,s em [0,6] com 5% em zero e resto log-uniforme em [0,01;6]; kappa,lambda em [0,20]; fisico).
2. Raios isolados: pontos com L_r<1 e min L_r. Massa isolada: pontos com alguma direcao natural passando nas razoes; distribuicao de q_max; q>=0,336 e q>=0,168.
3. Fronteira: min L_r entre q>=0,336, q>=0,168, q>=0,05 e q>0.
4. Refino Nelder-Mead (maxiter 300) com perda L_r+20*max(0,1-q/0,336) (40 se nenhuma direcao passa), a partir dos 20 melhores da tabela por essa perda mais 10 melhores por L_r puro (para provar raio isolado).
Viavel conjunto <=> L_r<1 e q>=0,336. Comparar com DYN12 (mesma fronteira: q=0,34 => L_r>=1,52).
## Antecipacao
Raios isolados: viavel em poucos pontos ou nenhum (com 8 sitios Sigma0 neutro (0,2,4) nao e equilatero; Sigma0 nao e identicamente zero). Fronteira massa-raio pior ou igual a DYN12; 0 viaveis conjuntos.
