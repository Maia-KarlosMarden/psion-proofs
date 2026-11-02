# P14 / WP-DYN: o hexagono dinamico - vertices que deslizam em passos de 60 graus e trocas antipodais (pre-registro)
Pergunta do autor: se os vertices se deslocam em angulos hexagonais e os psions trocam de posicao com o
antipoda, qual e a constituicao energetica?
## Fatos exatos usados (verificados antes)
(i) O antipoda k -> k+3 e o gerador Z2 da quiralidade (Z6 = Z3 x Z2): troca L <-> R.
(ii) Q_{k+3} - Q_k = +1 para k = 0,1,2 e -1 para k = 3,4,5 (exato): um salto de ROTULO ao antipoda muda a carga
     de +-1 (identidade carga-quiralidade Q = nR - 1 do paper).
(iii) A reacao CF-CC do paper, nu_mu(0,2,4) + n(0,1,5) -> e-(0,1,2) + p(0,4,5), e exatamente a troca dos
     psions antipodais 2 <-> 5 entre dois compostos (estrutura [T1] no paper).
## Modelo (rotulos mantem a carga; so a POSICAO e dinamica)
Estados: posicoes injetivas dos 3 rotulos da tripla nos 6 vertices (120 estados; exclusao de dois no mesmo vertice).
H = -eps sum_i cos(theta_{v_i} - theta_{k_i})        [vertice natural do rotulo k = 60k: o "relogio Z6"]
    - t  sum (passo de +-60 graus para vertice livre)   [deslizamento]
    - s  sum (salto de 180 graus para vertice livre + troca de posicoes de dois psions antipodais)  [antipoda]
Limite t = s = 0: hexagono estatico (controle: E0 = -3 eps, R_n = -0,154, R_Sigma- = 0,846).
Observaveis: E0, gap, peso do estado natural no fundamental, e R = <r^2>_t/<r^2>_p no estado fundamental
(mesma regra de WP-GEO: media de carga em torno do centroide, normalizada por Q_tot; sem normalizar se Q_tot = 0).
Varredura: t/eps em {0; 0,02; 0,05; 0,1; 0,2; 0,5; 1; 2; 5; 20}; s = 0 e s = t.
Limite delocalizado (eps -> 0): estado uniforme, R_n = 0 e R_Sigma- = 1 exatos (o resultado do bag).
## Antecipacao (por argumento, antes de rodar)
1) Custo de um passo de 60 graus: eps/2 por psion; de um salto antipodal: 2 eps (escala do relogio).
2) |R_n| diminui monotonicamente com t/eps, de 0,154 para 0; a dinamica nao aproxima de -0,164.
3) A troca antipodal quase nao altera o fundamental para s <~ eps; so importa para s >~ eps.
4) E0 desloca -O(t^2/eps) (perturbacao de 2a ordem). O estado natural pesa >90% so para t/eps <~ 0,1.
