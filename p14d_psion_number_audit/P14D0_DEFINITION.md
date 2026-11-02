# P14 / WP-D0: auditoria de cargas conservadas (pré-registro)

Base: seção 4 da especificação P14-D. Álgebra linear exata (Fractions), sem
parâmetros. Não usa η_B, ξ*, λ.

## Acréscimo declarado antes de rodar
O Teorema I.2 (J.6.2, [T1]) diz que o condensado de enrolamento ativa a ligação
|−1/2,L⟩ ↔ |+1/2,R⟩, isto é, um salto de um psion k = 2 ↔ k = 3, na fase
quebrada. Esse processo não estava na lista da seção 4 e é incluído aqui.

## Espaços
- Nível de constituintes: vetores de contagem de rótulos n ∈ Z⁶ (k = 0..5).
- Nível de compostos: contagens das 20 triplas de rótulos distintos (C(6,3)),
  observáveis (I = Σk mod 6 ∈ {0,3}, as 8 de WP-B) e não observáveis.
Uma carga conservada é um funcional linear c com c·Δ = 0 para todo processo.

## Processos (vetores Δ)
- G_UV(ψ̄ψ)²: não muda rótulos.
- JΣ (D1): ψ_k + ψ_l → ψ_{k+δ} + ψ_{l−δ}, k ∈ L = {0,1,2}, l ∈ R = {3,4,5},
  δ = ±1, com elemento de J± não nulo (o conjugado dá o reverso).
- JΣ (D1′): o mesmo, só componentes com ΔQ = 0, Q = (1,−5,−5,10,4,4)/9.
- Condensado: salto 2 ↔ 3 de um psion (fase quebrada).
- Rearranjo 2→2 entre triplas que preserva o multiconjunto de rótulos.
- Nível de compostos: JΣ intra-tripla (1→1) e entre triplas (2→2), salto do
  condensado dentro da tripla, e rearranjo; só alvos com rótulos distintos.
- Sphaleron: vetor desconhecido Δ_sph (P14-b), com ΔN_ψ ≠ 0 por hipótese
  (é a única fonte possível de N_ψ, porque a eq. (1) o conserva).

## Janelas
- S_sym: fase simétrica (fora da bolha), psions livres, JΣ ativo.
- S_brk: fase quebrada (dentro da bolha e depois do reaquecimento, já que
  T_RH < T_c em quase todas as linhas de WP-A), compostos com rearranjo, JΣ e
  condensado ativos (Γ/H ≫ 1 até T_f).
- Variante de sensibilidade S_brk′: sem o salto do condensado.

## Saídas
1. Bases racionais de K(S_sym) (constituintes, D1 e D1′) e de K(S_brk),
   K(S_brk′) (compostos).
2. Pertinência de N_ψ, M = Σ m_J, Q, n_R, B e L (esquemas S1–S3) a cada base.
3. Cargas candidatas: X ∈ K(S_brk) que não seja imagem de uma carga de K(S_sym)
   (isto é, violada na fase simétrica e conservada na quebrada). Só essas podem
   ser geradas na parede e sobreviver à lavagem.

## Critério (pré-registrado)
- VIÁVEL: existe candidata X. Reportar a condição sobre Δ_sph para X sobreviver
  (X·Δ_sph = 0) e se X tem projeção sobre B.
- INVIÁVEL: não existe candidata. Então, para qualquer Δ_sph, a assimetria
  gerada na parede de uma bolha Z6 é apagada no reaquecimento, e o transporte
  (WP-D) não precisa ser calculado.
- Q fora de K(S_brk) é registrado como achado independente (o condensado não
  conserva carga elétrica como simetria de Noether da eq. (1)).
