# P14 / WP-ABL: ablacao das premissas de massa (DYN12 e parafuso TWIST) (pre-registro)
DECLARACAO: pergunta do usuario ("e se a nossa presuncao de massa estiver errada?"). Registrado depois de TWIST. Adaptativo; as variantes sao escolhidas para ISOLAR cada premissa.
## Premissas de massa em uso (V0)
(P1) massa = funcional linear natural das esperancas <A>,<B>,<C> do estado fundamental, uma unica direcao (a,b,g) de mesmo sinal, comum a todos os hadrons; (P2) Lambda-n = 176,1 MeV > 0;
(P3) (Sigma0-Lambda)/(Lambda-n) = 0,437 +-30% (identifica Sigma0=(0,2,4)); (P4) (Sigma+-p)/(Lambda-n) = 1,425 +-30% (identifica Sigma+ como (1,3,5) ou (2,3,4), e liga o setor Sigma_k=9 ao 6);
(P5) escala natural: q = D/max|d_i| >= 0,336 (eps<=349 MeV, |w|<=1,5).
## Variantes (que premissas sao mantidas)
V0: P2,P3,P4,P5 (base). V1: sem P4 (so Sigma0). V2: sem P3 (so Sigma+). V3: sem P5 (escala livre; q>0). V4: so sinais (Sigma0>Lambda, Sigma+>p) + D>0 + P5. V5: so D>0 + P5. V6: sem massa (so raios).
## Procedimento
Para cada modelo (DYN12: 12 sitios fisico; TWIST: parafuso, 6 parametros) e cada variante: (1) tabela (DYN12: tabela de WP-SGN com lambda>0, 8000 pontos; TWIST: 30000 pontos LHS recomputados com as mesmas amostras de WP-TWIST);
 min L_r (perda de raios de EXH) entre pontos que satisfazem a variante; (2) refino Nelder-Mead (maxiter 400) dos 10 melhores por perda L_r+20*max(0,1-q/limiar) (40 se a variante nao e satisfeita; limiar 0,336 com P5, senao so q>0 => perda = L_r), mesmos limites de parametros que antes.
Viavel conjunto (variante) <=> L_r<1 e variante satisfeita. Reporta-se por variante o min L_r refinado e n de viaveis.
## Antecipacao
DYN12: raio puro viavel (V6: 0,79); V3 viavel (RIG, escala errada); V1, V2, V4, V5 NAO viaveis (L_r>1) porque o raio viavel exige t->0 e q ~ 0 => qualquer variante com P5 falha. TWIST: V6 inviavel (raio puro 1,15 > 1): nem a ablacao total da massa salva o parafuso.
Conclusao antecipada: o culpado no DYN12 e P5 (escala), nao P3/P4; no TWIST o problema e de raio (Lambda), nao de massa.
