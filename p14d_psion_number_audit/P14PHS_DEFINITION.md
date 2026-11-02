# P14 / WP-PHS: o anel Z6 como espaco de fase interno (pre-registro)
Leitura: os 6 rotulos sao os 6 angulos theta_k = 60 k graus da fase do condensado; o "raio de carga" e um
segundo momento de carga no espaco de fase, <r^2> = lambda^2 * sum_k Q_k d(theta_k, mu)^2 / Q_tot
(sem normalizar se Q_tot = 0), com UMA escala lambda. Nao existe mapa fase -> espaco derivavel da eq. (1)
(a interacao de contato G_UV nao depende de posicao); so os mapas abaixo, que o proprio paper distingue
(App. J.10.4: kappa vive na metrica plana da rede; G na metrica curva do circulo de fase; razao sinc(pi/6)=3/pi).
## Lista fechada
I1 metrica plana (corda, centroide euclidiano) = o hexagono de WP-GEO; CONTROLE: deve repetir R_n = -0,154, R_S- = 0,846.
I2 metrica curva (arco): d = distancia de arco, mu = media de Frechet dos tres angulos.
I3 corda com referencia na media circular: mu = arg(sum e^{i theta_k}), d = |e^{i theta_k} - e^{i mu}|.
Triplas equilateras (1,3,5), (0,2,4): sum e^{i theta} = 0 e a media de Frechet nao e unica -> I2 e I3 indefinidas
(registrado como "sem previsao", nao como falha).
## Alvos (os de WP-GEO8)
n: -0,164 a 10%; Sigma-: 0,86 +- 0,21; Sigma0: |R| <= 0,15; Lambda: (0; 0,25]; Sigma+: [0,80; 1,15].
## Antecipacao (feita a mao antes de rodar - nao cega)
Em arcos consecutivos I2 e I3 dao <r^2> proporcional a (Q_tot - Q_centro): R_n = -0,200 (22% do medido),
R_S- = 0,80, Lambda = +1,0, (2,3,4) = -0,2. Logo I2/I3 nao melhoram o hexagono; I1 continua melhor no nucleon.
