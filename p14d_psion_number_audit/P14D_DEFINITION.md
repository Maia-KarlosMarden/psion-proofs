# P14-D: pressão como segunda variável na bariogênese da PH (pré-registro)

Registrado antes de rodar. Entradas proibidas: η_B observado (só na comparação
final), ξ* = 434, Ω_DM.

## Achados de leitura (antes de calcular)

1. K.10: as paredes Z6 interpolam a fase do condensado entre vácuos vizinhos com
   rigidez f_π² e gap m_Z6, ambos de escala MeV (f_π = 0,4994 MeV,
   m_Z6 = 2,996 MeV, σ^{1/3} = 0,423 MeV). Com barreira ~ m_Z6² f_π² ≪ T⁴ para
   T ≫ MeV, as paredes só existem como defeitos para T ≲ MeV.
2. K.10 Layer 2: viés κ_CP = (m_e/Λmeta)², o mesmo da bariogênese; aniquilação
   em t_ann ~ σ/ε ≈ 4 anos. O texto diz "a assimetria que faz a matéria mata as
   paredes": no paper a pressão do viés é consequência, não causa.
3. Z6 é global (C-I8.5): há paredes, não vórtices. Um vórtice de rotação tem
   pressão central menor; pressão central maior exige tensão de parede (Laplace
   2σ/R) ou pressão de vácuo (domínio de falso vácuo comprimindo plasma).
4. Na transição super-resfriada de P14-c, a pressão que move as bolhas é
   ΔV = λv⁴/8, e α = ΔV/ρ_rad(T_nuc) = (T_eq/T_nuc)⁴.

## Pré-checagem de escalas (este registro)

Para cada época candidata, calcular:
- Γ_rearr/H com Γ_rearr = c_Γ g⁴ T⁵/Λmeta⁴ (taxa de contato quatro-férmions,
  G_UV = g²/Λmeta², g² = 2π²/3; c_Γ ∈ [10⁻³, 1] cobre fatores de espaço de fase
  não derivados) e H de radiação (g* = 106,75 acima de 100 GeV; 10,75 em MeV;
  3,36 abaixo de 0,5 MeV).
- Atividade do sphaleron (T_max = 1,2×10¹¹ GeV de WP-C; inativo abaixo de ~130 GeV).
- Épocas: formação das paredes Z6 (T = m_Z6), aniquilação (t = 4 anos, T pela
  relação tempo–temperatura de radiação), e as nove linhas resolvidas de WP-A
  (T_nuc, T_RH = T_eq).
- Para as linhas de WP-A: α = (T_eq/T_nuc)⁴; diluição D_TI = (T_nuc/T_RH)³; e a
  cota de bariogênese espontânea em regime lento
      η_cota = (45 κ /(2π² g*)) · ε_ψ · (π/3) · (Γ_rearr/T)|_{T_nuc} · D_TI,
  κ = 1, para c_Γ = 10⁻³ e 1. Isto é uma COTA DE VIABILIDADE com fatores O(1)
  livres, não o cálculo de transporte.

## Critério da pré-checagem

- Versão "paredes Z6": viável só se Γ_rearr/H > 10⁻³ ou sphaleron ativo em
  alguma época em que as paredes existem. Caso contrário, excluída.
- Versão "bolhas da transição P14-c": segue para o cálculo de transporte
  (WP-D, a especificar) se η_cota com c_Γ = 1 alcançar ≥ 10⁻¹⁰ em pelo menos
  uma linha; caso contrário, registrar que a pressão de vácuo, ao produzir
  inflação térmica, reduz a assimetria abaixo do observado em todo o intervalo.

## Transparência

Ao redigir, fiz de cabeça uma estimativa preliminar (η ~ 10⁻¹⁷ para a linha
λ = 1/25) e a registro aqui antes de rodar, para que o resultado da pré-checagem
possa ser comparado com ela.
