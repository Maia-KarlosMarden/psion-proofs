# P14 / WP-ROB: massa sem Sigma0 + viabilidade robusta (DYN12) (pre-registro)
DECLARACAO: segue WP-ABL (P3, a razao (Sigma0-Lambda)/(Lambda-n), nao e sustentada; V2/V4/V5 chegaram a 0,99 mas na borda e frageis: 0/60 sob +-10%). Registrado depois de ABL. Adaptativo.
## Criterio de massa (V2 de WP-ABL; Sigma0 fora, hipotese: desdobramento hiperfino de spin que um modelo sem spin nao gera)
Direcao natural (mesmo sinal) de pesos sobre <A>,<B>,<C> (4000 direcoes Dirichlet +/-): D=M_L-M_n>0; (M_Sp-M_p)/D = 1,425 +-30% (Sp = (1,3,5) ou (2,3,4)); escala q=D/max|d_i|>=0,336. SEM razao nem sinal para Sigma0.
## Raios (inalterados): n -0,164 +-0,0164; Sigma0 |R|<=0,15 (raio de Sigma0 continua exigido: so a MASSA do Sigma0 e retirada); Lambda [0,05;0,25]; Sigma+ [0,80;1,15].
## Perda nominal L_c(z) = L_r + 20*max(0, 1 - q/0,336) (40 se nenhuma direcao passa). Viavel nominal <=> L_c<1.
## Viabilidade ROBUSTA
Conjunto fixo de 6 perturbacoes (semente 40): y = z + 0,10*max(|z_i|, 0,05)*u, u_i ~ U[-1,1] por coordenada (z em (t,s,kappa,lambda)), recortado ao dominio fisico [0,6]x[0,6]x[0,20]x[0,20].
Perda robusta L_rob(z) = max(L_c(z), max_j L_c(y_j)). Viavel robusto <=> L_rob<1.
## Procedimento
Dominio: DYN12 (12 sitios), lambda,kappa>=0. Partidas (11): as 8 melhores da tabela de WP-SGN (lambda>0) por L_c nominal; os otimos de V2 e V5 de WP-ABL; o ponto robusto de EXH (0,001; 0,069; 0,26; 0,11).
Nelder-Mead sobre L_rob (maxiter 200). Validacao independente: 100 perturbacoes novas (semente 41, mesmas regras) em cada otimo; reporta-se a fracao com L_c<1.
Sucesso <=> L_rob<1 e fracao de validacao >= 90%.
## Antecipacao
Nenhuma partida chega a L_rob<1 (melhor entre 1,2 e 1,8); fracao de validacao < 50%. Se houver sucesso, teste de acaso e verificacao de escala dos pesos antes de qualquer afirmacao.
