# P14 / WP-C2: a fase do vácuo é física sem contagem dupla do bond? (pré-registro)

Motivo: em WP-SGN o Hamiltoniano das triplas tinha um bond real |2⟩↔|3⟩ (h₁ de WP-SIM) somado
ao termo do condensado J⁺Σ⁺ + J⁻Σ⁻ com fase. Mas o J.6.2 diz que o bond é GERADO por esse mesmo
condensado. Contagem dupla suspeita.

## Parte 1 (algébrica, exata)
H_k = diag|m| + g O_JS + λ (e^{iθ_k} X + e^{−iθ_k} X†), X = (J⁺P_L)(J⁺P_R), sem bond separado,
g = 1, λ = 1/6. Verificar no espaço de Fock (64) e na base de triplas (20):
(a) [M̂, O_JS] = 0, [M̂, diag|m|] = 0, [M̂, X] = 2X (M̂ = Σ m_k n_k);
(b) U = e^{−iθ_k M̂/2}: ‖U H_k U† − H_0‖ = 0 para k = 0..5 (fase removível por U(1)_M).

## Parte 2 (dinâmica diagnóstica)
Acrescentar uma quebra explícita de U(1)_M de força δ (bond real ΔM = 1 entre k = 2 e 3, como em
h₁): H_k(δ) = H_k + δ·bond. Rede quântica de WP-LOOP (L2, ω = 10, γ = 1, c_Γ = 1e-3; usa psions
livres, só como detector de CP-ímpar), k = 1, δ ∈ {0, 0,01, 0,03, 0,1, 0,3, 0,5}.
Saída: Δp(δ); expoente de escala d ln|Δp|/d ln δ.

## Parte 3 (tamanho na PH)
A única quebra explícita de U(1)_M registrada no paper é a separação dos vácuos de K.10,
κ_CP = (m_e/Λ_meta)². Estimar ε ~ κ_CP^p com o expoente p da Parte 2 e η ~ (η/ε)_max·ε com
(η/ε)_max = 0,16 (WP-COMP).

## Controles
C1 hermiticidade. C2 δ = 0: |Δp| no piso (< 1e-15). C3 k = 0: Δp = 0 para todo δ.

## Antecipação registrada
Parte 1: fase exatamente removível (‖·‖ < 1e-12) em todos os vácuos. Parte 2: Δp = 0 em δ = 0 e
Δp ∝ δ² (o invariante de CP mais simples é Im(t² V*)). Parte 3: ε ~ κ_CP² ~ 1e-62 (ou ~1e-31 se o
expoente for 1); η ≲ 1e-32: a fonte de CP do vácuo não sustenta a bariogênese na PH atual. O
resultado positivo de WP-SGN em PB/PC é atribuído à contagem dupla.
