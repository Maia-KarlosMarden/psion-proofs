# P14 / WP-J1: a assimetria de J.1 nas duas conjugações do paper (pré-registro)

Ontologia do paper (A.5, J.1): matéria e antimatéria são tuplas de psions; a
antipartícula de uma tupla é a sua imagem por uma conjugação de índices.
Duas conjugações aparecem no texto e são testadas em paralelo:
- C_A: k → (6 − k) mod 6 (A.5, L.4.2).
- C_B: k → (k + 3) mod 6, o Z2 ⊂ Z6 que troca L ↔ R (J.1: LLL ↔ RRR).
Sem η_B, ξ*, λ.

## Ingredientes (do paper)
- Hamiltoniano de um psion de J.6.3: h = diag(|m_k|) + t(|2⟩⟨3| + |3⟩⟨2|),
  m_k = k − 5/2, t = 1/2 (Teorema I.3). Unidades: as de h.
- Vértice J·Σ na convenção D1 (geradores completos, A_a = J_a P_L, B_a = J_a P_R,
  com o conjugado), como operador de dois corpos no espaço de Fock de 6 modos:
  O_JS = ½ Σ_a [(c†A_a c)(c†B_a c) + h.c.]. Variante D1′ (ΔQ = 0) também.
- Intensidade relativa g do vértice em relação a h: não fixada pelo paper em
  unidades de h; reportar g ∈ {0,1; 1} (sensibilidade, 2 valores declarados).
- Setor de 3 psions: 20 triplas (base de Slater).

## Passo 2: paridade
Para cada conjugação C (matriz de permutação P), decompor h e O_JS em partes par
e ímpar: X± = (X ± PXPᵀ)/2. Reportar ‖X−‖/‖X‖. J.1 afirma O_JS ímpar sob Z2.

## Passo 3: diferença de energia entre tupla e conjugada
H = H₀ + g O_JS no setor N = 3. Para cada tripla observável t e C(t):
(a) ΔE₀ = ⟨t|H₀|t⟩ − ⟨C(t)|H₀|C(t)⟩ (só h);
(b) ΔE₁ = ⟨t|H|t⟩ − ⟨C(t)|H|C(t)⟩ (primeira ordem em g);
(c) energias dos autoestados exatos de H com maior sobreposição com t e C(t).
Reportar também Q(t) e Q(C(t)) com Q = (1,−5,−5,10,4,4)/9.

## Critério (pré-registrado)
- J.1 confirmado na conjugação C se: (i) O_JS tem parte ímpar não nula sob C;
  (ii) o próton (0,4,5) e o elétron (0,1,2) têm energia menor que as suas
  conjugadas por C; (iii) Q(C(t)) = −Q(t) para os dois. O critério (iii) é o de
  conjugação de carga; se falhar, a conjugação não descreve antimatéria
  observável, qualquer que seja (i)–(ii).
- Reportar o resultado de cada conjugação, sem escolher.
