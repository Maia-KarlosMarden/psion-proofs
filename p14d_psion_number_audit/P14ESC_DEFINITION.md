# P14 / WP-ESC: restricao de escala absoluta no DYN12 (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: motivado pelo adendo de WP-RIG (D=0,016 => eps ~ 10,7 GeV, incoerente com eps ~ 163-349 MeV). Registrado depois de RIG. Nao cego.
## Restricao
Massa M_h = eps * (a<A>_h + b<B>_h + g<C>_h) em MeV, com eps (energia de relogio) em [163, 349] MeV (anel: 163 MeV; correcao relativistica: ~349 MeV) e pesos naturais (mesmo sinal)
com max|w_i| <= 1,5 (w=(1,1,1) e a energia simples). Razoes de massa: como em RIG (D=M_L-M_n>0, (S0-L)/D em 0,437 +-30%, (Sp-p)/D em 1,425 +-30%; n-p fixado por constante c).
Escala: Lambda-n = 176,1 MeV => D(w)=176,1/eps em [0,505; 1,080]. Para uma direcao d (soma 1), w=s*d com s=(176,1/eps)/D_d; exige s*max|d_i|<=1,5. Condicao mais folgada (eps=349): q=D_d/max|d_i| >= 0,336.
q e invariante de escala da direcao: q e a "alcance de massa" de um ponto dinamico.
## Procedimento
1. Tabela de WP-SGN (lambda>0, 8000 pontos), 4000 direcoes (Dirichlet +/-): para cada ponto, q_max = max q entre direcoes que passam nas razoes; reporta distribuicao de q_max e quantos pontos tem q_max>=0,336.
2. Mesmo teste com limite alternativo max|w|<=3 (q>=0,168) como sensibilidade.
3. Refino Nelder-Mead (maxiter 250), perda = L_r + 20*max(0, 1 - q/0,336) se alguma direcao passa nas razoes; 40 se nenhuma. Partidas: 20 melhores da tabela por perda e o otimo de EXH.
Viavel conjunto <=> L_r<1 e q>=0,336. Dominio fisico t,s em [0,6]; kappa,lambda em [0,20].
## Antecipacao
q_max da tabela << 0,336 em quase todos os pontos (0 com q>=0,336 na tabela; ou < 1%); apos refino 0 viaveis conjuntos. Se houver pontos com q >= 0,336, eles terao raios fora da janela.
