# P14 / WP-ESC2: varredura fina ao redor do quase-acerto de escala natural (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: dirigida pelo melhor ponto de WP-ESC (t,s,k,l)=(0,899; 0,014; 0; 0,79), L_r=1,55, q=0,341. Registrado depois de ESC. Nao cego.
## Grade
t em {0,5; 0,6; 0,7; 0,8; 0,9; 1,0; 1,1; 1,2; 1,3; 1,4; 1,5; 1,75; 2,0} (13); s em {0; 0,01; 0,02; 0,04; 0,08; 0,12; 0,2} (7); kappa em {0; 0,05; 0,1; 0,2; 0,4; 0,8} (6);
lambda em {0,4; 0,5; ...; 1,4} (11 passos de 0,1): 13*7*6*11 = 6006 pontos. Modelo DYN12, estado fundamental.
## Criterios (de WP-ESC, inalterados)
Radios: L_r<1 (perda de EXH). Escala/massa: existe direcao natural (4000 direcoes Dirichlet +/-) que passa nas razoes (D>0; (S0-L)/D 0,437 +-30%; (Sp-p)/D 1,425 +-30%)
com q = D/max|d_i| >= 0,336 (eps<=349 MeV, |w|<=1,5). Sensibilidade: q>=0,168 (|w|<=3).
Etapas: (1) grade: contar L_r<1 e L_r<1,2 entre q>=0,336 e entre q>=0,168. (2) Nelder-Mead (maxiter 300) a partir dos 10 melhores pontos da grade com perda L_r+20*max(0,1-q/0,336).
Viavel conjunto <=> L_r<1 e q>=0,336. Reportar tambem o ponto de menor L_r com q>=0,168.
## Antecipacao
L_r minimo com q>=0,336 fica em 1,2-1,6 (n e Lambda nao entram juntos na janela com a escala); 0 viaveis conjuntos; com q>=0,168 pode haver L_r<1.
