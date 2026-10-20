# P14 / WP-OPS: termos de CP simples feitos de operadores da PH (pré-registro)

Pergunta: algum termo de violação de CP com 1 ou 2 parâmetros, construído só com
operadores que a PH já tem, reproduz as 24 restrições da calibração (WP-CAL)?
Se sim, a bariogênese pela montagem vira previsão (só ε ajustado). Se nenhum se
sair melhor que direções sorteadas, a receita de WP-CAL é acomodação.

## Rede e restrições
As de WP-CAL: c_Γ = 1e-3, lei L2 (receita) e, como sensibilidade, L1. Matriz resposta
G (26 × 69) recalculada como em WP-CAL. Restrições A1–A4 (Δp = Δe = D, Δν = −2D,
16 tuplas exóticas e 5 psions livres em zero): vetor b (24).

## Candidatos (cada um é uma direção; o termo da Lagrangiana não é escrito, só o
## padrão de CP que ele induz nas taxas)
Operadores diagonais na base de tuplas: Q (carga, soma de Q9), M (soma de m_k),
N_R (número de índices ≥ 3), Z (soma dos índices mod 6), E (diagonal de h + J·Σ).
Operadores não diagonais: h e O = J·Σ (matrizes 20×20 reais de WP-NOW).
- Transições (antissimétricas, como exige CPT): K = [A, B] para A ∈ {Q, M, N_R, Z, E}
  e B ∈ {h, O}, mais K = [h, O]. 11 candidatos. Padrão y_e = K_ij nas arestas,
  projetado nas circulações (núcleo de B_E).
- Rearranjos ab → cd: c_r = f(a)f(b) − f(c)f(d) com f ∈ {Q, M, N_R, Z, E} (o
  produto não é aditivo; somas aditivas se cancelam porque o rearranjo conserva
  os psions). 5 candidatos, projetados no núcleo da incidência por pares.
Total: 16. A parte de cada candidato removida pela projeção (proibida pela
unitariedade) é descartada e sua fração é reportada.

## Medida
Ajuste por mínimos quadrados de ε (1 termo) ou (ε₁, ε₂) (2 termos, 120 pares)
para G u ≈ b; resíduo relativo ρ = ‖G v − b‖/‖b‖ e ε_ef da solução.
Linha de base: 2000 direções gaussianas sorteadas no espaço de 69 dimensões
(semente 0), para 1 termo e 2000 pares para 2 termos.

## Critério
Sucesso (previsão possível): algum candidato de 1 ou 2 termos com ρ < 0,1 e
ε_ef < 1. Estrutura (sem sucesso): melhor candidato com ρ abaixo do 1º percentil
da linha de base. Caso contrário: nada especial.

## Antecipação registrada
Espero ρ ~0,9–1 para todos (um ou dois números não zeram 16 tuplas exóticas) e
nenhum candidato fora da linha de base. Ou seja, espero que a receita de WP-CAL
não venha de um termo simples.
