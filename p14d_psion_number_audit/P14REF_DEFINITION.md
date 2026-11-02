# P14 / WP-REF: refino local conjunto raio+massa (DYN12, massa reponderada) (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: segue WP-SGN (23% dos pontos passam em massa com pesos livres, mas L de raios >= 3,28 na grade grossa). Nao cego.
## Perda conjunta
Parametros z=(t,s,kappa,lambda,a,b,g); dominio t,s em [0,6]; kappa em [0,20]; lambda em [0,20] (fisico); pesos conforme a variante.
Raios: L_r = max(|R_n+0,164|/0,0164, |R_S0|/0,15, |R_L-0,15|/0,10, min_j|R_Spj-0,975|/0,175) (perda de EXH).
Massa: M_h = a<A>_h + b<B>_h + g<C>_h; D=M_L-M_n; v1=|(M_S0-M_L)/D / 0,437 - 1|/0,3; v2=min_j|(M_Spj-M_p)/D / 1,425 - 1|/0,3; se D<=0: L_m=5+|D|; senao L_m=max(v1,v2).
(c comum {n,S0,L} vs {p,Sp} fixa n-p=1,29 MeV e nao entra na perda.) Perda: L_c = max(L_r, L_m). Viavel <=> L_c<1.
## Variantes
A: pesos livres em [-3,3]^3 (vacuo conhecido: controle de WP-SGN: alvos aleatorios de massa passam em 54% dos pontos). B: pesos nao-negativos e naturais em [0,1,5]^3 (sem inversao de sinal; (1,1,1) e a energia simples).
## Procedimento
Partidas: 12 melhores pontos por L_c da tabela de WP-SGN (lambda>0), para cada variante, mais o otimo de EXH (t,s,k,l)=(0,001; 0,069; 0,26; 0,11) com pesos (1,1,1).
Nelder-Mead (maxiter 250). Reporta-se min L_c por variante, e para L_c<1: R, pesos, D, razoes e se |a,b,g| em [0,1,5] com mesmo sinal (natural).
## Regra de interpretacao (fixada antes)
Sucesso de A sozinho NAO e evidencia (vacuo). Sucesso de B (pesos naturais, sem sinal invertido) seria o primeiro candidato raio+massa e exigiria teste de acaso proprio antes de qualquer afirmacao.
## Antecipacao
A: min L_c entre 0,5 e 1,0 (pesos livres resolvem a massa e o raio EXH ja viavel); B: min L_c > 1 (nenhum ponto viavel).
