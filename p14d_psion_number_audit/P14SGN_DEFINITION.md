# P14 / WP-SGN: (T1) sinal de J.Sigma invertido e (T2) massa como funcional reponderado (DYN12) (pre-registro)
DECLARACAO: pedido do usuario ("Teste os") apos WP-MASS. Registrado depois de conhecer MASS e do diagnostico (E_L<E_n, E_n<E_p, E_Sp2<E_p em todos os pontos). Adaptativo; nao cego.
Leitura previa do v86: eq.(1) tem -(psibar J psi)(psibar Sigma psi) (atrativo), mas nenhuma derivacao fixa o sinal da ligacao hadronica. Portanto lambda<0 e EXPLORATORIO (nao derivado).
## Tabela
16000 pontos LHS: t,s em [0,6] (5% em zero, resto log-uniforme [0,01;6]); kappa em [0,20] (fisico); lambda em [-20,20]. Modelo DYN12 (12 sitios, 1320 colocacoes, fundamental).
Para cada hadron: R2, E_gs, e as esperancas <A> (relogio), <B> (Coulomb), <C> (J.Sigma) no fundamental (E_gs = <A>+kappa<B>+lambda<C> + hopping).
## T1: lambda<0 (kappa>=0)
Criterios de WP-MASS: D=E_L-E_n>0; R1=(E_S0-E_L)/D em 0,437 +-30%; R2=(E_Sp-E_p)/D em 1,425 +-30% (Sp1 ou Sp2); c fixado por n-p=+1,29 MeV. Depois menor L de raios (perda de EXH; viavel <=> L<1).
## T2: massa reponderada (lambda>0, kappa>=0 fisico)
Massa M_h = a <A> + b <B> + g <C> + c_h, com pesos (a,b,g) em {-3,-2,-1,-0.5,0,0.5,1,2,3}^3 (729) livres e universais (os mesmos para todos os hadrons) e c comum a {n,S0,L} vs {p,Sp}. Funcao de onda e raios vem da dinamica (t,s,kappa,lambda).
Criterios: mesmos D'>0, R1', R2' com M no lugar de E. Reporta-se por ponto se existe ALGUM peso que passa; depois, entre esses, menor L de raios.
Controle de vacuidade (obrigatorio): 3 pesos livres vs 3 razoes (D', R1', R2') => quase sempre soluvel. Mede-se a fracao de pontos com algum peso passando para (i) alvo real e (ii) 200 alvos de massa aleatorios
(r1 alvo U[0,1;1,0], r2 alvo U[0,5;3,0], mesmas janelas relativas de 30%). Se (i) ~ (ii), o T2 nao e informativo.
## Antecipacao
T1: D>0 passa em muitos pontos (inversao do sinal), mas R1&R2&raios nao passam juntos (0 com L<1); T2: muitos pontos com algum peso passando (vacuo), e entre eles os raios continuam L>1; (i) ~ (ii).
