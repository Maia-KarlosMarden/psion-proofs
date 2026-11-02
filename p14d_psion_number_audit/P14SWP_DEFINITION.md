# P14 / WP-SWP: varredura abrangente do dodecagono dinamico (pre-registro)
DECLARACAO DE ADAPTATIVIDADE: motivada pelos 4 pontos individuais de Lambda vistos em WP-DYN12, que nao coincidiram com n e Sigma+.
Mesmo modelo de WP-DYN12 (12 sitios, 1320 colocacoes, saltos +-1/+-2 sitios com -t, troca/salto antipodal com -s, fundamental via eigsh).
## Grade (fina)
t, s em {0, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3.5} (12 valores cada); kappa, lambda em {0, 0.5, 1, 2, 3, 5, 8, 12} (8 cada):
12*12*8*8 = 9216 pontos, 5 hadrons (+p).
## Criterios (inalterados)
Raios: n -0,164 +-0,0164; Sigma0 |R|<=0,15; Lambda [0,05;0,25]; Sigma+ [0,80;1,15] em (1,3,5) ou (2,3,4).
Massas (apenas entre os viaveis em raio): E_n-E_p>0; E_L-E_p>0, escala u=177,4/(E_L-E_p) MeV; Sigma+ - p = 251 +-30%; Sigma0 - Lambda = 77 +-30%
(energias do estado fundamental). Teste de acaso na MESMA grade densa: 2000 alvos aleatorios (centros U[-0,6;2,2]).
## Regra de interpretacao (fixada antes)
Um ponto viavel so conta como evidencia se P_acaso da grade <= 0,05. Se P_acaso > 0,05, a viabilidade nao distingue o modelo de ajuste generico.
Pontos em ambigua degenerescencia (t=s=0) excluidos (NaN).
## Antecipacao
Alguns pontos viaveis em raio podem aparecer (grade 41x mais densa), mas P_acaso sobe; massas falham (o sinal E_n>E_p e E_L>E_p nao e garantido).
