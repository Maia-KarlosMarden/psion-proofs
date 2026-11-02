# P14 / WP-TWIST: parafuso Z3 x Z2 (prisma torcido) dinamico (pre-registro)
DECLARACAO: registrado depois do fracasso de WP-PRISM (produto direto Z3 x Z2, todas as triplas com uma classe Z3 de cada => parte angular degenerada). Hipotese do usuario: Z2 em estrutura MISTA com Z3.
Adaptativo; nao cego quanto ao fracasso do prisma e da fronteira massa-raio dos poligonos.
## Modelo (familia de dois parametros geometricos)
6 sitios (camada s//3, classe Z3 s%3). Camada L (k=0,1,2): angulos 0,120,240 graus, z=0. Camada R (k=3,4,5): angulos 120*(k%3)+delta, z=h. Rotulo k tem sitio natural k.
 delta=0: prisma (WP-PRISM); delta=60 graus: antiprisma; delta=180 graus com h->0: estrela de dois triangulos invertidos (hexagono com L = triangulo equilatero, antipodas a 180 graus).
 A inversao de quiralidade (Z2) leva a classe j da camada L para a classe j da camada R COM torcao delta e translacao h: Z2 e Z3 nao comutam geometricamente (parafuso), ao contrario do prisma.
Dinamica, termos e criterios: identicos a WP-PRISM (hopping -t no triangulo da mesma camada; -s inversao de quiralidade na mesma classe e troca de antipodas; fundamental por diagonalizacao densa 120x120;
 A = -sum [cos(phi_sitio - phi_natural) + (1 se mesma camada, -1 se outra)]/2; B Coulomb 3D; C J.Sigma 3D).
Parametros livres: t, s in [0,6]; kappa, lambda in [0,20]; h in [0,4;3] log-uniforme; delta uniforme em [0, 180 graus].
## Criterios (inalterados)
Raios L_r<1 (n -0,164 +-0,0164; S0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15]); massa: direcao natural (mesmo sinal), D=M_L-M_n>0, (S0-L)/D 0,437 +-30%, (Sp-p)/D 1,425 +-30%, q=D/max|d_i|>=0,336 (sensibilidade 0,168).
## Procedimento
1. Tabela LHS de 30000 pontos (6 parametros). 2. Raios isolados, massa isolada (q_max), fronteira min L_r para q>=0,336, 0,168, 0,05, >0; distribuicao em delta dos melhores. 3. Refino Nelder-Mead (maxiter 400, 6 parametros): 20 partidas conjuntas (perda L_r+20*max(0,1-q/0,336), 40 se sem massa) e 10 por L_r puro.
Viavel conjunto <=> L_r<1 e q>=0,336. Em caso de sucesso, teste de acaso proprio antes de qualquer afirmacao.
## Antecipacao
delta=0 reproduz o prisma (q indefinido em todos); para delta>0 aparecem pontos com razoes de massa; fronteira massa-raio entre a do DYN12 e a do octogono (min L_r em q>=0,336 entre 1,5 e 5); 0 viaveis conjuntos (probabilidade de sucesso conjunto ~15%).
