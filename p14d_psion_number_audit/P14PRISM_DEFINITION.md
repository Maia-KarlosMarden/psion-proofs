# P14 / WP-PRISM: prisma triangular Z3 x Z2 dinamico (pre-registro)
DECLARACAO: registrado depois do fechamento do DYN12 e de WP-OCT. Geometria escolhida por ARGUMENTO DE SIMETRIA, nao por ajuste: o PH tem Z6 = Z3 x Z2, antipoda k->k+3 = gerador da quiralidade (Z2),
cargas Q_k = q(k mod 3) + floor(k/3) com q=(1,-5,-5)/9 (App. J.6.4: Q_{k+3}-Q_k = +-1). Logo a posicao Z3 e angular (0,120,240 graus) e a quiralidade e a camada. Isso distingue o prisma do hexagono (que mistura as duas
em um angulo unico). Ainda e adaptativo: conhece-se a tensao massa-raio dos poligonos.
## Modelo
6 sitios: sitio s = (camada s//3, angulo 120 s%3 graus), coordenadas (cos, sin, h*camada); R do triangulo = 1; h = distancia entre camadas, PARAMETRO LIVRE (log-uniforme em [0,4; 3]). Sitio natural do rotulo k = sitio k (L=camada 0, R=camada 1).
120 colocacoes ordenadas por tripla. H = diag + hopping:
 diag = A + kappa*B + lambda*C com A = -sum_a [cos(dphi_a) + (1 se mesma camada, -1 se outra)]/2 (relogio, eps=1), B = sum q q/d (3D), C = -sum_{pares L-R} 1/d (3D).
 hopping: -t salto no triangulo (+-120 graus, mesma camada, sitio livre); -s inversao de quiralidade (mesmo angulo, outra camada, sitio livre) e troca de dois psions antipodais (mesmo angulo, camadas diferentes).
 Estado fundamental (diagonalizacao densa 120x120). Raios: <r^2> carga-ponderado em torno do centroide / Q_tot (Q_tot=0: sem normalizar) em 3D, razoes a <r^2>_p.
## Criterios (inalterados de WP-ESC/OCT)
Raios: L_r<1 (n -0,164 +-0,0164; S0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] em (1,3,5) ou (2,3,4)).
Massa: M_h = a<A>+b<B>+g<C>, direcao natural (mesmo sinal), 4000 direcoes Dirichlet +/-: D=M_L-M_n>0, (S0-L)/D em 0,437 +-30%, (Sp-p)/D em 1,425 +-30%; escala natural q=D/max|d_i|>=0,336 (eps<=349 MeV, |w|<=1,5); sensibilidade q>=0,168.
## Procedimento
1. Tabela LHS 30000 pontos (t,s log-uniformes em [0,01;6] com 5% em zero; kappa,lambda em [0,20]; h em [0,4;3]).
2. Raios isolados; massa isolada (q_max); fronteira min L_r para q>=0,336, 0,168, 0,05, >0.
3. Refino Nelder-Mead (maxiter 400, 5 parametros) com perda L_r+20*max(0,1-q/0,336) (40 se nenhuma direcao passa): 20 melhores da tabela; e 10 melhores por L_r puro.
Viavel conjunto <=> L_r<1 e q>=0,336.
## Antecipacao
Raios isolados viaveis apos refino; fronteira massa-raio melhor que o octogono mas ainda com L_r>1 em q>=0,336 (0 viaveis conjuntos; probabilidade de sucesso conjunto ~20%).
Em um sucesso, o resultado exigira teste de acaso proprio (alvos aleatorios) antes de qualquer afirmacao.
