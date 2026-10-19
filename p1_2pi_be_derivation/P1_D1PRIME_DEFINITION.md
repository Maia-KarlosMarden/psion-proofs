# P1 — D1′: vértice e condensado que conservam carga elétrica (pré-registro)

Entrada nova no registro. Não apaga D1 (vértice completo J·Σ, peso violador de
carga 0,598) nem D3/D3′. Gatilho: o teste de carga do passo 4 disparou a
cláusula "viola conservação de carga" da D1. O motivo é físico e independe do
alvo 32/25.

## 1. Definição

Cargas (as mesmas do passo 4): Q = (1, −5, −5, 10, 4, 4) para k = 0..5, j = 5/2.

a) **Vértice.** Cada operador A_a = J_a P_L e B_a = J_a P_R é decomposto em graus
   de carga: A_a^(q) contém só os elementos |k′⟩⟨k| com Q[k′] − Q[k] = q. O
   vértice D1′ é V′ = Σ_a Σ_q A_a^(q) ⊗ B_a^(−q), isto é, só entram componentes
   com ΔQ total nulo. A parte J_z⊗J_z é diagonal e entra inteira.
b) **Condensado.** s é Hermitiana suportada nos estados livres (k ≥ 1) e só tem
   elementos fora da diagonal entre estados de mesma carga: pares permitidos
   (1,2) e (4,5). O condensado é neutro.
c) Resto como em D3′: W[s] = (Tr s)² − Tr(s²)/N_D + σ[Re Σ Tr(A^(q) s)Tr(B^(−q) s)
   − Re Σ Tr(A^(q) s B^(−q) s)/N_D], λ_max no emparelhamento de Frobenius,
   W_crit′ = n_livres·λ_max, Z′ = 2^n_livres / W_crit′, b_e′ = b0/(Z′·N_c).
d) Quatro ramos (σ, N_D) reportados, nenhum adotado.

## 2. Limitações declaradas antes de rodar

- As cargas Q só estão definidas para j = 5/2. **Não há teste de família em j**
  para D1′, então o critério forte (forma fechada em j) não se aplica a ela.
- Esquema: a razão-limite W[s*]/Tr(s*²) → λ_max é verificada de novo com os dois
  reguladores em j = 5/2.

## 3. Critério de resultado e regra de parada (fixados antes de rodar)

- "Compatível": |Z′/(32/25) − 1| ≤ 1% em algum ramo. Mesmo assim, como a família
  em j não pode ser testada, isso é compatibilidade numérica, não evidência;
  com 12 tentativas acumuladas (D3 4 + D3′ 4 + D1′ 4) a chance de acerto
  casual em ±1% fica da ordem de 10%. Nenhum ramo é adotado.
- Qualquer outro resultado é registrado como achado negativo.
- Regra de parada: após D1′ resta no máximo uma redefinição (D3″) antes de
  concluir, em relatório, que be não sai do 2PI nas escolhas feitas.

## 4. Transparência

Os valores com condensado restrito e vértice completo (D1, D3-N) já eram
conhecidos do passo 4 (W = 18,41 e 18,14 para σ = +1). Os valores de D1′ com λ_max
não foram calculados antes deste registro.
