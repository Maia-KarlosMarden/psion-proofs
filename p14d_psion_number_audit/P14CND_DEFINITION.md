# P14 / WP-CND: linha de base condicionada nos raios para as previsoes P-A e P-B (pre-registro)
DECLARACAO: segue WP-OOS (p0 = 0/82881 NAO condicionado nos raios; O3 passa em P-A e P-B). Caixa de amostragem escolhida DEPOIS de ver os otimos de HYB3 (adaptativo). Pergunta: dentro da regiao que ja satisfaz os raios e o V0, com que frequencia as previsoes P-A e P-B passam por acaso?
## Procedimento
Caixa 1: t em [0; 0,5], s em [0; 0,8], kappa em [0; 0,5], lambda em [0; 0,6]; Caixa 2 (sensibilidade): t em [0; 1,0], s em [0; 1,6], kappa em [0; 1,0], lambda em [0; 1,2]. Amostragem uniforme (LHS) de 16000 pontos por caixa.
Modelo, regra (dinamico <=> contem 2 ou 3; p, n estaticos), raios e V0 identicos a HYB3/OOS: filtrar pontos com L_r<1 (n -0,164 +-0,0164; S0 |R|<=0,15; L [0,05;0,25]; Sp [0,80;1,15]) e com ao menos 1 direcao natural aprovada em V0 (q>=0,336, razoes a 30%).
Para os pontos que sobrevivem, avaliar P-A (R_Sm em [0,65;1,07] e massa rel. n em [181;335] MeV) e P-B (massa de (3,4,5) rel. n em [205;381] MeV) em todas as direcoes aprovadas (4000 direcoes Dirichlet +/-).
Metricas: (i) n de pontos viaveis V0; (ii) fracao dos pontos viaveis com ALGUMA direcao que passa P-A&P-B (p_ponto); (iii) fracao de pares (ponto, direcao) V0 que passa P-A&P-B (p_par); idem so P-A e so P-B.
Interpretacao (fixada antes): se p_ponto >= 0,20, o desempenho de O3 e comum dentro da regiao viavel (nao e evidencia); se p_ponto <= 0,05, O3 e raro (evidencia fraca, ainda sujeita a look-elsewhere).
## Antecipacao
p_ponto entre 0,05 e 0,25 (P-B facil: massas largas; P-A limitada pelo raio do Sigma-, quase sempre perto de 1,0-1,4); O3 nao e raro.
