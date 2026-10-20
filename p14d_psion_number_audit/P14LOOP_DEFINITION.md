# P14 / WP-LOOP: uma fase em J·Σ, dinâmica quântica da montagem (pré-registro)

Pergunta: um termo concreto de CP, de um só parâmetro, seleciona sozinho p, e⁻ e ν̄
(as restrições de WP-CAL), sem padrão heurístico? A unitariedade sai da própria
mecânica quântica, não de projeção.

## Regra adicionada (fora da eq. (1), declarada)
J·Σ → Re(e^{−iφ} O) = cos φ · O_JS + sin φ · O_A, com O = Σ_a (ψ̄J^aP_Lψ)(ψ̄J^aP_Rψ)
na base de Fock de WP-J1, O_JS = (O + O†)/2 (o da eq. (1)) e O_A = (O − O†)/2i.
Na base de tuplas O_A é imaginário puro e hermitiano (8 entradas), e
max |Im H_ab H_bc H_ca| = 3,8 em φ = 0,3: há fase de CP invariante.
Hamiltoniano das tuplas H(φ) = h + g Re(e^{−iφ}O), g = 1; antimatéria: H(φ)*.

## Dinâmica (substitui as transições clássicas W de WP-NOW)
Por setor, matriz densidade ρ (20×20) das tuplas e abundâncias dos 6 psions:
dρ/d(ln) = r(x) { −iω[H, ρ] − γ(ρ − diag ρ) + diag(k_b P_aP_bP_c)
  − ½{diag(k_b/K), ρ} − ½{diag(p̄), ρ} + diag(Y_eq,t²) + rearranjos }
- r(x) = Γ/H de WP-NOW; ω e γ (frequência coerente e decoerência por colisões, em
  unidades da taxa de contato) são escalas que a PH não fixa: grade declarada.
- Ligação ψψψ → |t⟩ (estado de base), dissociação, aniquilação t t̄ (p̄ = população
  da antitupla correspondente; no setor de antimatéria, o inverso) e rearranjos
  clássicos sobre as populações (só na diagonal), como em WP-NOW; psions como em
  WP-NOW (ligação e aniquilação ψψ̄).
- Equilíbrio (ρ = Y_eq,t·I, psions em 1/12) é estacionário para qualquer φ.

## Grade
lei ∈ {L2, L1} × c_Γ = 1e-3 × ω ∈ {0,1, 1, 10} × γ ∈ {0,1, 1, 10}, φ = 1e-3, δ = 0.
18 corridas. Sensibilidade no melhor (ω, γ) de cada lei: c_Γ ∈ {1e-5, 1}.

## Medida
Vetor das 24 restrições de WP-CAL: d = (Δp, Δe, Δν₁+Δν₂, 16 exóticas, 5 psions
livres). Resíduo com um só fator de escala s (o tamanho de sen φ):
ρ = min_s ‖s d − b‖/‖b‖, b = (1, 1, −2, 0…). ε necessário: o φ que dá η = 6,1e-10.
Comparar com WP-OPS (1 termo, direções sorteadas): 1º percentil 0,926 (L2) e 0,787 (L1).

## Controles
C1 φ = 0: assimetria < 1e-12. C2 |N_ψ| < 1e-12. C3 Δ(−φ) = −Δ(φ) (1e-6 relativo).
C4 do equilíbrio com φ = 0,3 e taxas constantes: assimetria < 1e-12.
Linearidade: Δ/φ entre φ = 1e-3 e 1e-2 (reportar).

## Critério
Sucesso: ρ < 0,1. Estrutura: ρ abaixo do 1º percentil da linha de base. Senão: nada.

## Antecipação registrada
Espero ρ entre 0,8 e 1 (sem estrutura): a fase em J·Σ não conhece a regra de níveis
e não deve silenciar as 16 tuplas exóticas. Espero que a assimetria exija 3+ tuplas
acopladas (análogo de Jarlskog) e se concentre na componente grande do grafo.
