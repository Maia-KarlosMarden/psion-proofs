# P14 / WP-D0′: cogênese por conservação de N_ψ (pré-registro)

Registrado antes de extrair massas e números de corpos dos 39 estados escuros.
Entradas proibidas: η_B, ξ*, λ, Ω_DM/Ω_b observado (só na comparação final).

## Premissa estrutural (de D0 e da retificação de 19/10)
A eq. (1) conserva N_ψ exatamente, e nada feito de psions o viola. Logo a
assimetria só pode ser uma SEPARAÇÃO de N_ψ entre setores que deixam de trocar
psions: visível (compostos de 3 corpos), neutrinos, escuro (N ≥ 9).
Com processos entre setores ativos, a única carga conservada é N_ψ total; depois
do congelamento deles, N_ψ de cada setor é conservado separadamente. Isto é
estrutural e não depende de taxas.

## Hipóteses declaradas
- H1: N_ψ total = 0 (sem condição inicial).
- H2: o setor de neutrinos tem N_ψ líquido zero. (Alternativa não testada: os
  neutrinos carregam parte do balanço; isso torna a relação abaixo uma cota.)
- H3: a matéria escura relíquia é feita de antipsions (estados conjugados dos 39;
  mesma massa por CPT).
- H4: massas e estabilidade dos 39 estados como no repositório
  (`P_yield_ratio.py`, passo D.6: lei de cascata com b_e = 325/1152, A = 13/12,
  p = 5/6, Λ = 1,2662e12 GeV; estabilidade por `_stable`). Isto herda o
  estatuto de postulado de b_e (P1) e a nota "em revisão" do Teorema KS.

## Relação (sem parâmetro livre)
    N_ψ(escuro) = −N_ψ(visível)  ⇒  n_DM · N_DM = ν_vis · n_b
    Ω_DM/Ω_b = ν_vis · (m_DM / N_DM) / m_N
    ν_vis = 3 + 3 (n_e/n_b),  n_e/n_b = 1 − Y_p/2,  Y_p = 0,245  ⇒  ν_vis = 5,6325
    m_N = 0,9389 GeV (média próton-nêutron)

## Regra de escolha do relíquia (fixada antes de ver as massas)
- R1 (principal): dentro do setor escuro, a T ≪ massas, processos que conservam
  N_ψ concentram o número no estado de menor massa por psion, min(m_i/N_i).
- R2 (sensibilidade): densidade igual em todos os 39 estados:
  m/N efetivo = Σm_i / ΣN_i.
- R3 (sensibilidade): o estado mais leve, min(m_i), se diferente de R1.
Sensibilidade a Λ: repetir R1 com Λ = 1,2609e12 (valor do paper; m ∝ Λ).

## Critério
Compatível se R1 ficar a ±10% de Ω_DM/Ω_b observado (Planck 2018,
Ω_c h² / Ω_b h² = 0,1200/0,02237 = 5,364). Tentativas: 3 regras (R1, R2, R3).
Compatível em R1 conta como relação estrutural não ajustada, condicional a
H1–H4; não fixa η_B nem o mecanismo de separação (CP e não equilíbrio), que é o
passo seguinte.

## Transparência
Ao redigir, fiz de cabeça a conta inversa: o observado exige m/N ≈ 0,89 GeV, isto
é, m_DM ≈ 8 GeV para N = 9 e ≈ 13 GeV para N = 15. Sei também, de uma saída
anterior do repositório, que há "duas espécies degeneradas a 14,3 GeV", sem saber
o N delas. Registro isso antes de rodar.
