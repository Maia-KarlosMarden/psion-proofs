# P14 / WP-HYB3: HYB2 sem Sigma- no laco, refino robusto e teste de acaso (pre-registro)
DECLARACAO: registrado depois de HYB2 (refino nao concluido; ponto nominal de V0 em (0,303; 0,195; 0,044; 0,193), 85% estavel, Sigma- falha). Adaptativo.
## Modelo (regra de HYB2)
Estaticos: p (0,4,5), n (0,1,5). Dinamicos (DYN12, 12 sitios, fundamental): Sigma0 (0,2,4), Lambda (1,2,3), Sigma+ (1,3,5) e (2,3,4). Sigma- (0,1,2) FORA do laco e dos criterios (tripla do eletron, inerte a J.Sigma); sera reportado a parte
no otimo (como dinamico e como estatico) sem entrar em nenhum criterio. Isto e uma decisao tomada depois de ver que Sigma- falha e trava o solver.
## Criterios (V0 de WP-HYB, inalterados)
Raios: n -0,164 +-0,0164; Sigma0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] (opcao (1,3,5) ou (2,3,4)). Massa: direcao natural (mesmo sinal), 4000 direcoes, D=M_L-M_n>0, (M_S0-M_L)/D 0,437 +-30%,
(M_Sp-M_p)/D 1,425 +-30%, q=D/max|d_i|>=0,336. Perda L_c = L_r + 20*max(0,1-q/0,336). Robusta L_rob = max(nominal, 6 perturbacoes +-10% semente 40); validacao 100 perturbacoes semente 41 (regras de WP-ROB).
## Procedimento
1. Tabela LHS 8000 pontos (semente 43, igual a HYB2) com R e esperancas <A>,<B>,<C> por hadron (salvas).
2. Teste de acaso no nivel de tabela (mesmo procedimento para alvo real e aleatorios): 300 conjuntos de alvos aleatorios: centros de raio U[-0,6;2,2] (n, S0, L, Sp), alvos de massa r1 U[0,1;1,0], r2 U[0,5;3,0], mesmas janelas relativas. P_acaso(a) = fracao com algum ponto da tabela com L_c<1.
   P_acaso(b): igual, mas com o centro do raio de n fixado no valor real (-0,164): o R_n estatico e uma previsao sem parametros e nao deve ser contado como acaso.
   Alvo real na mesma tabela: n de pontos com L_c<1.
3. Refino Nelder-Mead (maxiter 100) sobre L_rob a partir das 10 melhores da tabela; validacao.
Sucesso <=> L_rob<1 e validacao >=90% E P_acaso(b) <= 0,05.
## Antecipacao
Refino acha ao menos 1 ponto robusto (prob. ~50%); P_acaso(a) < 0,02 (n estatico fixa R_n); P_acaso(b) entre 0,03 e 0,2 (criterio de massa com 3 pesos livres e frouxo). Sigma-: raio > 1,1 e massa muito acima de Sigma0 nos otimos.
