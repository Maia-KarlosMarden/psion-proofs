# P14: mecanismo de bariogênese da PH (B−L, sphalerons, ξ*)

Três pacotes de trabalho pré-registrados (definição com SHA-256 antes de rodar;
`P14_*_REGISTRY.txt` guarda hashes e carimbos UTC, inclusive das correções de código).

| Pacote | Pergunta | Definição | Script | Saída |
| --- | --- | --- | --- | --- |
| WP-B | De onde vem B−L ≠ 0? | `P14_WPB_DEFINITION.md` | `P14_wpb.py` | `P14_wpb_output.txt` |
| WP-C | O sphaleron vale em T_nuc? Ccore é reprodutível? | `P14_WPC_DEFINITION.md` | `P14_wpc.py` | `P14_wpc_output.txt` |
| WP-A 1 | Qual transição permite ξ* ≫ 1? | `P14_WPA1_DEFINITION.md` | `P14_wpa1.py` | `P14_wpa1_output.txt` |
| WP-A 2 | Bounce no potencial de Coleman–Weinberg (P14-c) | `P14_WPA2_DEFINITION.md` | `P14_wpa2.py` | `P14_wpa2_output.txt` |

Resultados, em resumo:
- A dinâmica Z6 conserva só N_ψ; recombinações 2→2 violam B, L e B−L. B−L ≠ 0 sai
  da partição térmica entre famílias, com sinal dependente da identificação dos
  estados neutros LLR (P14-a). O fator que substitui ε_conf = 4/9 é c_s(B−L),
  3 a 34 vezes menor.
- ξ* = 434 era calculado de η_B (`xi_star_from_etaB`); a eq. (38) do v85 é identidade.
  Ccore (eq. 14) não é reprodutível. O 8/23 de §13.3 é o fator de sphaleron para
  dois dubletos de Higgs. A eq. (33) não tem a diluição de entropia (fator 27,3).
- Sphalerons do SM: ativos em T_nuc ≈ 3×10⁹ GeV (Γ/T³H ≈ 54), fora de equilíbrio
  em Λmeta. Escala de emergência do SU(2)_L não especificada no paper (P14-b).
- NJL em campo médio (σ₀ = 3 MeV) é de 2ª ordem com T_c = 9,6 MeV; Potts Z3 dá
  ξ* ≤ 1,08. Só um potencial conforme (P14-c) permite ξ* ≫ 1.
- Com P14-c e cinco λ pré-registrados, ξ* fica entre 3,6 e 39 (máximo em
  λ = 1/(2j)²), T_c ≈ 10¹¹ GeV, e há inflação térmica em todos os casos
  resolvidos. A eq. (33) não descreve esse regime.

Pré-requisitos: Python 3, numpy, scipy. `P14_wpa2.py` leva ~50 min.
