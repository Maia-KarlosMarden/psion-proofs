# P14 / WP-CP: auditoria de violação de CP na PH (pré-registro)

Pergunta: a PH tem uma fonte de CP que possa gerar assimetria líquida a partir
de um estado inicial CP-simétrico? Cálculo exato; sem η_B, ξ*, λ.

## Convenções
Bilinear escalar com matriz interna T: sob C, ψ̄Tψ → ψ̄Tᵀψ; sob P, inalterado
(J, Σ comutam com γ^μ, K.9); logo sob CP, T → Tᵀ = T* (T hermitiana).
A eq. (1) é CP-invariante se o seu tensor interno de quatro férmions é real na
base k (a menos de uma rotação unitária interna).

## Testes
A1. Montar o tensor interno de (ψ̄ψ)² e de (ψ̄J^aψ)(ψ̄Σ^aψ) nas convenções D1
    (geradores completos, D1-i com o conjugado) e D1′ (graus de carga com ΔQ=0).
    Medir max|Im|. Real ⇒ nenhuma violação explícita de CP na eq. (1).
A2. ε_ψ: com hopping complexo t = |t|e^{iφ} na ligação |−1/2,L⟩↔|+1/2,R⟩,
    calcular ε_ψ(φ) a partir do Hamiltoniano de J.6.3 (diag(|m|) + t|L⟩⟨R| + h.c.)
    e verificar ε_ψ(φ) = ε_ψ(−φ). Se par em φ, ε_ψ é CP-par: mede a força da
    mistura L/R, não uma assimetria entre psion e antipsion.
A3. Condensado de enrolamento w J⁺Σ⁺ + w* J⁻Σ⁻: sob CP, w → w*. Nos vácuos Z6,
    w_k = |w| e^{2πik/6}. Listar Im(w_k^n) para n = 1, 2, 3: vácuos k = 0, 3 são
    CP-invariantes; os demais violam CP espontaneamente, com sinais opostos em k e
    6−k (vácuos CP-conjugados).
A4. Viés entre vácuos (K.10, Layer 2): ΔC(s) = C(0) − Σ_k m_J(k) m_J(k+s).
    Verificar ΔC(s) = ΔC(6−s). Se simétrico, vácuos CP-conjugados são degenerados
    pelo próprio viés do paper; para qualquer população p_k = f(E_k),
    Σ_k p_k Im(w_k^n) = 0 exatamente.
A5. κ_CP = (m_e/Λmeta)²: número real positivo; registrar valor e que, por si, não
    é uma fase.

## Critério (pré-registrado)
- SEM FONTE LÍQUIDA DE CP: A1 real, A2 par, A4 simétrico. Então, com estado
  inicial CP-simétrico, a PH (eq. 1 + vácuos e viés de K.10) só produz assimetrias
  de sinais opostos em domínios CP-conjugados, com soma zero. A assimetria
  observada seria condição inicial, como N_ψ.
- FONTE EXPLÍCITA: A1 com parte imaginária; reportar o tamanho.
- FONTE ESPONTÂNEA LÍQUIDA: A4 assimétrico; reportar Σ p_k Im(w_k^n).
Rotas não testadas aqui (só listadas): condição inicial CP-assimétrica; termos
complexos fora da eq. (1) (por exemplo um termo θ); seleção de vácuo por
mecanismo não contido em K.10.
