# P14 / WP-HYB: hibrido estatico/dinamico (pre-registro)
DECLARACAO: segunda parte do pedido ("faca um depois o outro"). Registrado depois de ROB. Adaptativo; a regra abaixo e fixada ANTES de rodar e nao e escolhida hadron a hadron.
## Regra de selecao (principio, nao por estado)
Um hadron e DINAMICO <=> sua tripla contem o par ligado (2,3) (a unica ligacao real |2>-|3> do App. J.6.3, h1[2,3]=0,5): Lambda=(1,2,3) e Sigma+=(2,3,4). Os demais (p=(0,4,5), n=(0,1,5), Sigma0=(0,2,4)) sao ESTATICOS: psions fixos nos vertices naturais do anel (sitio 2k do dodecagono).
Para a Sigma+ so a opcao (2,3,4) vale (a (1,3,5) nao contem o par).
## Modelo
Dinamicos: DYN12 (12 sitios, 1320 colocacoes, relogio eps=1 + kappa*Coulomb + lambda*J.Sigma + saltos -t (+-1,+-2 sitios) e -s (troca/salto antipodal), estado fundamental). Parametros (t,s,kappa,lambda) em [0,6]x[0,6]x[0,20]x[0,20].
Estaticos: colocacao natural, esperancas <A>=-3, <B>, <C> da colocacao natural; raios R_n=-0,154 e R_S0=0 (identidade do arco, independentes dos parametros), R relativo ao <r^2>_p estatico.
## Criterios
Raios: n -0,164 +-0,0164; Sigma0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] (perda L_r de EXH).
Massa (duas versoes): V0: D=M_L-M_n>0, (M_S0-M_L)/D 0,437 +-30%, (M_Sp-M_p)/D 1,425 +-30%, q=D/max|d_i|>=0,336; V2: sem a razao do Sigma0 (como WP-ROB). Direcao natural (mesmo sinal), 4000 direcoes.
Perda nominal L_c = L_r + 20*max(0,1-q/0,336) (40 se nenhuma direcao passa). Robusto: max sobre nominal e 6 perturbacoes de +-10% (semente 40, regras de WP-ROB); validacao em 100 perturbacoes novas (semente 41).
## Procedimento
(1) Tabela LHS 8000 pontos (t,s,kappa,lambda): raios (so Lambda e Sigma+ variam), q, perda por versao. (2) Nelder-Mead (maxiter 200) sobre L_rob a partir das 10 melhores de cada versao. (3) Validacao.
Sucesso <=> L_rob<1 e validacao >=90%.
## Antecipacao
Raios: viaveis (n e Sigma0 estaticos passam por construcao; Lambda e Sigma+ alcancam as janelas no regime deslocalizado). V2: robusto viavel (criterio frouxo). V0 (com Sigma0): inviavel; (Sigma0 estatico e n estatico tem energia natural comparavel e Lambda dinamico fica mais ligado: D<0 ou q baixo). Probabilidade de sucesso de V0 ~25%.
