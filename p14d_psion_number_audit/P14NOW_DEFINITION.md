# P14 / WP-NOW: a PH como está, simulada (pré-registro)

Ontologia O− (WP-ONT): antipartícula = antipsions, N_ψ = −3, carga −Q.
Dinâmica: só a eq. (1) (h de J.6.3 e J·Σ D1), CP de campo exata (WP-CP).
η_B observado só na comparação final.

## Etapa A (1,38e11 → 1e8 GeV; processos de contato, Γ/H como em WP-SIM)
Espécies: 6 psions, 6 antipsions, 20 tuplas, 20 antituplas (52).
Processos: ligação ψψψ ↔ t e ψ̄ψ̄ψ̄ ↔ t̄ (pesos L1/L2 de WP-SIM, K = 10);
transições internas t ↔ t′ e t̄ ↔ t̄′ com |⟨t′|h + g e^{±iφ} O_JS|t⟩|² (sinal +
para matéria, − para antimatéria: CP conjuga a fase); rearranjo 2→2 em cada setor;
aniquilação ψ_k ψ̄_k ↔ radiação e t t̄ ↔ radiação (peso 1, equilíbrio relativístico).
Condição inicial: Y_ψk = (1+δ)/12, Y_ψ̄k = (1−δ)/12.
Grade: c_Γ ∈ {1e-5, 1e-3, 1} × lei ∈ {L1, L2} × δ ∈ {0, 1e-8} × φ ∈ {0, π/4};
g = 1. 24 corridas.
Saídas: assimetria líquida de prótons Δ_p = Y(0,4,5) − Y((0,4,5)‾), de elétrons
Δ_e = Y(0,1,2) − Y((0,1,2)‾), e N_ψ líquido.

## Etapa B (congelamento da aniquilação p p̄; T de m_p até m_p/1000)
Boltzmann padrão: Y_p − Y_p̄ = Y_B constante (de A),
dY_p̄/dx = −(λ/x²)[Y_p̄(Y_p̄ + Y_B) − Y_eq²], x = m/T,
Y_eq = 0,145 (g/g*) x^{3/2} e^{−x}, λ = 0,264 (g*/√g*) M_Pl m ⟨σv⟩,
g = 2, g* = 17,25, ⟨σv⟩ = 1/m_π² (taxa hadrônica medida, entrada externa),
m = 0,938 GeV (sensibilidade: 0,610 GeV, a escala do nível 8).
Saída: n_p/n_γ e n_p̄/n_γ finais (n/n_γ = 7,04 Y).

## Controles
C1 δ = 0, φ = 0: Δ_p = Δ_e = 0 e N_ψ líquido = 0 a precisão de máquina.
C2 N_ψ líquido conservado em toda corrida (desvio < 1e-12 relativo a δ).
C3 φ ≠ 0 com δ = 0: assimetria zero (taxas de árvore não distinguem e^{±iφ}
   com amplitudes reais; é a exigência de fase forte).

## O que se espera (antecipação registrada)
δ = 0: η = 0 e uma relíquia simétrica de prótons e antiprótons ~10⁻¹⁹–10⁻¹⁸ de
n_γ cada. δ ≠ 0: η final proporcional a δ (a PH só redistribui o que for posto
no início). φ não muda nada.
