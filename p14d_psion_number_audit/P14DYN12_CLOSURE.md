# P14 / Fechamento do DYN12 (dodecagono dinamico) -- resumo registrado
Registrado em 2026-10-31 apos WP-ESC2. Todos os WPs abaixo foram pre-registrados (SHA-256 + UTC) antes de rodar; as antecipacoes erradas estao em cada _REGISTRY.txt.

## Modelo
12 sitios, 1320 colocacoes ordenadas por tripla, H = relogio(eps=1) + kappa*Coulomb + lambda*J.Sigma + saltos (-t: +-30/60 graus; -s: troca/salto antipodal), estado fundamental.
Alvos de raio (em unidades de <r^2>_p): n -0,164 +-0,0164; Sigma0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15]. Massa: razoes (Sigma0-Lambda)/(Lambda-n)=0,437 e (Sigma+-p)/(Lambda-n)=1,425 a 30%, n-p por constante.

## Linha do tempo (WP: resultado)
- DYN12: 0/225 viaveis; Lambda alcanca a janela em 4 pontos, nunca junto com n e Sigma+.
- SWP (9216 pts): 0 viaveis; n+Lambda juntos em 4; P_acaso 0,003.
- EXH: ponto viavel EM RAIO (L=0,796; t=0,001 s=0,069 k=0,26 l=0,11), estavel +-10%; massas com sinal errado (m_n<m_p, m_L<m_p).
- RND: P_acaso continuo 0/80, mas pipeline truncado nao recupera o proprio alvo real (limite enviesado).
- MFIRST: kappa>=0 da m_n<m_p e m_L<m_p em 12000/12000 pontos (B_p=+0,275, B_n=-0,155 => E_n-E_p=-0,43 kappa). Sinal fixado pela atribuicao de cargas, nao por principio de exclusao/incerteza.
- MASS: constantes por hadron (mu*Sum(k-5/2)+nu*Q) nao alteram E_L-E_n; D>0 em 0/12000 (J.Sigma domina: corr(E,lambda)~-0,9).
- SGN: lambda<0 nao resolve (E_n>E_p 0/8000); massa reponderada com pesos livres passa em 23% mas e quase vacua (alvos aleatorios 54%).
- REF: refino com penalidade colapsa (artefato); sem conclusao.
- RIG: 4/21 viaveis conjuntos (raio + razoes de massa com pesos naturais), mas cone de pesos estreito (2/4000 direcoes), D=0,016 => eps~10,7 GeV (30-65x fora de 163-349 MeV).
- ESC: com eps 163-349 MeV e |w|<=1,5, 129/8000 pontos tem alcance de massa q>=0,336; refino 0/21 viaveis conjuntos (melhor L_r=1,55).
- ESC2: q>=0,336 => min L_r=1,52-1,60; q>=0,168 => min L_r=1,078; 0 viaveis conjuntos.

## Resultado central: fronteira massa-raio (exclusao observada)
| alcance de massa q | melhor L de raios | origem |
| 0,016 | 0,79 (viavel) | RIG; escala 30-65x errada |
| 0,243 | >= 1,078 | ESC2 (|w|<=3) |
| 0,340 | >= 1,52 | ESC2 (eps<=349 MeV, |w|<=1,5) |
Monotona nos pontos testados: quanto melhor a escala de massa, pior o raio. Raio viavel exige localizacao (t->0); separacao de massa na escala certa exige t~0,9 (delocalizacao).
Hipotese NAO provada (a testar): e uma relacao tipo incerteza entre localizacao (raio) e separacao de energia por mistura cinetica (massa). Nao ha demonstracao analitica nem varredura global.

## Status e limites
- DYN12 FECHADO como mecanismo conjunto raio+massa no dominio fisico [0,6]^2 x [0,20]^2 (kappa,lambda>=0) com busca em grade, hipercubo latino, evolucao diferencial e Nelder-Mead locais. Nao e prova de ausencia global.
- Raios isolados: reproduzidos (L=0,796) com massas erradas. Massas isoladas: reproduziveis com pesos livres (vacuo). Conjunto em escala natural: nao encontrado.
- Nada disto altera o status de A1/A2 (baryogenese) nem das conclusoes do v86 alem de L.4.7; fica como resultado negativo para o paper se o usuario decidir inclui-lo.
- Pendente: octogono (8 sitios, antipoda = +4; rotulos k=0,1,2,3,4,5 em sitios 0,1,2,4,5,6; sitios 3 e 7 vagos) como proxima geometria.
