# P14 / WP-RIG: massa como restricao rigida com pesos naturais; minimizar raios (DYN12) (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: corrige o desenho de WP-REF (plateau da penalidade). Registrado depois de REF. Nao cego.
## Natural
Massa M_h = a<A>_h + b<B>_h + g<C>_h. Como as razoes de massa sao invariantes de escala, so a DIRECAO (a,b,g) importa. Direcao NATURAL: todos os componentes com o mesmo sinal
(octante positivo ou negativo, sem inversao relativa entre relogio, Coulomb e J.Sigma). Teste: 4000 direcoes amostradas (Dirichlet(1,1,1), 2000 positivas e 2000 negativas); um ponto e "massa-natural"
se alguma direcao amostrada passa em D=M_L-M_n>0, R1=(M_S0-M_L)/D em 0,437 +-30%, R2=(M_Sp-M_p)/D em 1,425 +-30% (Sp1 ou Sp2). (n-p fixado pela constante c.)
## Etapas
1. Tabela de WP-SGN (lambda>0; 8000 pontos): contar pontos massa-naturais; menor L_r (perda de raios de EXH) entre eles.
2. Controle de vacuidade: 40 alvos de massa aleatorios (r1 alvo U[0,1;1,0], r2 alvo U[0,5;3,0]) -> fracao de pontos massa-naturais (real vs aleatorio).
3. Refino: Nelder-Mead (maxiter 250) minimizando L_r com filtro rigido (perda = L_r se massa-natural, senao 20), partindo dos 20 melhores massa-naturais por L_r e do otimo de EXH. Dominio fisico: t,s em [0,6]; kappa,lambda em [0,20].
Viavel conjunto <=> L_r<1 e massa-natural. Reporta-se a direcao natural encontrada (a,b,g)/norma e a similaridade com (1,1,1) (energia simples).
## Regra de interpretacao
Se o alvo real for MAIS facil que a media dos aleatorios, a restricao natural nao e informativa. Um viavel conjunto exigiria teste de acaso proprio antes de qualquer afirmacao.
## Antecipacao
Massa-naturais: ~5-15% dos pontos (octante unico reduz o vacuo de 23%); min L_r entre eles > 1 apos refino (0 viaveis conjuntos).
