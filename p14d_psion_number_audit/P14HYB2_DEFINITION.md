# P14 / WP-HYB2: hibrido com regra 'contem 2 ou 3' (pre-registro)
DECLARACAO: motivado pelo adendo de WP-HYB (massa ~ <A>: deslocalizacao => mais pesado; mas Sigma0 estatico quebra o triploto Sigma). Registrado depois de HYB. Adaptativo; a regra e fixada antes e e por principio:
## Regra
DINAMICO <=> a tripla contem o psion 2 ou o 3 (os dois psions do par ligado J.6.3, h1[2,3]): Lambda (1,2,3), Sigma0 (0,2,4), Sigma- (0,1,2), Sigma+ em ambas as opcoes (1,3,5) e (2,3,4). ESTATICOS: p (0,4,5) e n (0,1,5), as unicas triplas sem 2 nem 3.
Modelo, parametros, raios, massa V0/V2, perda nominal e robusta, validacao: identicos a WP-HYB. Sigma- (0,1,2) e simulado (inerte a J.Sigma: C=0) e reportado: raio R_Sm (obs 0,86 +-0,21) e massa.
Sigma0 e agora dinamico: |R_S0|<=0,15 deixa de ser automatico (agora e um criterio a atender).
Criterios adicionais REPORTADOS (nao na perda): desdobramento do triploto Sigma: |M_Sm-M_S0|/D e |M_Sp-M_S0|/D (obs 0,03 e 0,02) na direcao aprovada; ordem Sigma0>Lambda.
## Procedimento
Tabela LHS 8000 pontos; NM sobre L_rob a partir das 10 melhores de V0 e de V2; validacao 100 pertubacoes novas (semente 41). Sucesso <=> L_rob<1 e validacao >=90%.
## Antecipacao
V2 robusto viavel (>=1 de 10 partidas); V0 (com razao Sigma0) viavel com prob. ~30%; |R_S0|<=0,15 atendido; Sigma- radius em [0,65;1,07] em minoria dos pontos; triploto Sigma quase degenerado (Sigma0, Sigma-, Sigma+ dinamicos com <A> proximos).
