# P14 / WP-LJS: derivar a escala de ligacao barionica Lambda_JS do termo J.Sigma (pre-registro)
Alvo: Lambda_JS = hbar c / a = 0,163 GeV (a = 1,212 fm, fixado pelo proton no hexagono; nota DT).
Equivale a ln(Lambda_meta/Lambda_JS) = 29,68, ou degrau l* = 8,37 da cascata.
Criterio de sucesso: previsao sem parametro livre com |ln(previsto/alvo)| <= 0,69 (fator 2).
## Parte A: derivacao dinamica em 4D (gap de NJL no canal J.Sigma)
Termo de quatro fermions com corte Lambda_meta, canal atrativo com autovalor lambda_ch de J.Sigma
(p, n: -10). Equacao de gap de campo medio (corte euclidiano 4D):
 1 = (G_eff N Lambda^2 / 4 pi^2) [1 - (M^2/Lambda^2) ln(1 + Lambda^2/M^2)].
Saida independente da normalizacao de G: a sintonia Delta = 1 - G_c/G_eff exigida para M/Lambda = alvo.
Se Delta << 1e-10, a escala nao e derivada, e sintonizada (hierarquia por potencia, nao exponencial).
## Parte B: transmutacao exponencial (so possivel se J.Sigma for marginal), lista FECHADA com
constantes do paper (b0 = 13/6, b_e = 325/1152, N_c = 6, g^2 = 2 pi^2/3, S = 212,677):
 B1 um laco 4D com o b0 da PH: Lambda_meta exp(-8 pi^2/(b0 g^2)).
 B2 duas insercoes do condensado: Lambda_meta exp(-S/g^2).
 B3 degrau barionico inteiro do P14 (l = 8): Lambda_meta exp(-8/b_e).
 B4 tipo Kondo com o casimir do spin interno j = 5/2 (rho J = C g^2/8 pi^2, C = j(j+1)):
    Lambda_meta exp(-8 pi^2/(C g^2)).
Look-elsewhere: 4 formas; com expoentes espalhados em ~0-100, a chance de uma cair a +-0,69 por
acaso e ~4 x 1,4/100 ~ 5%.
Antecipacao: A exige sintonia ~1e-24 (nao deriva). B: nenhuma acerta; B2 erra por ~14, B3 por ~3,7.
