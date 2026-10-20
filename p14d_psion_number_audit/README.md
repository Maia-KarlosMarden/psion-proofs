# P14-D: pressão, lavagem e conservação do número de psions

Pacote pré-registrado (definição com SHA-256 antes de rodar; `*_REGISTRY.txt`
guarda hashes e carimbos UTC). Complementa `p14_baryogenesis/` (PR #10).

| Etapa | Pergunta | Definição | Script | Saída |
| --- | --- | --- | --- | --- |
| Pré-checagem | Em que época a pressão pode gerar assimetria? | `P14D_DEFINITION.md` | `P14D_scales.py` | `P14D_scales_output.txt` |
| Adendo pós-hoc | A assimetria sobrevive ao reaquecimento? | — | `P14D_washout_addendum.py` | `P14D_washout_addendum_output.txt` |
| WP-D0 | Que cargas cada fase conserva? | `P14D0_DEFINITION.md` | `P14D0_audit.py` | `P14D0_audit_output.txt` |
| WP-D0′ | A conservação de N_ψ fixa Ω_DM/Ω_b? | `P14D0b_DEFINITION.md` | `P14D0b_cogenesis.py` | `P14D0b_output.txt` |
| WP-CP | A eq. (1) viola CP de campo? | `P14CP_DEFINITION.md` | `P14CP_audit.py` | `P14CP_output.txt` |
| WP-J1 | J.1 nas duas conjugações do paper | `P14J1_DEFINITION.md` | `P14J1_test.py` | `P14J1_output.txt` |
| WP-M | Tupla e conjugada têm a mesma massa no mesmo nível? | `P14MASS_DEFINITION.md` | `P14MASS_test.py` | `P14MASS_output.txt` |
| WP-SIM | Rede de Boltzmann de formação de tuplas | `P14SIM_DEFINITION.md` | `P14SIM_network.py` | `P14SIM_output.txt` |

Resultados:
- Paredes Z6 (K.10) existem de T ≈ 3 MeV até ~4 anos; nada viola B nessa época. Excluídas.
- Bolhas da transição P14-c: cota de viabilidade passa, mas cobre 13 ordens de
  grandeza e escala como T_nuc⁷; o rearranjo está em equilíbrio no reaquecimento
  (Γ/H = 32 a 1,4×10⁵) e lava B, L e B−L.
- A eq. (1) conserva o número de psions exatamente. Nada feito de psions o viola,
  inclusive um "sphaleron" emergente. Logo B+L = N_ψ/3 é conservado e só B−L
  muda: o inverso do Modelo Padrão. A conversão por sphaleron (28/79) usada em
  WP-B/WP-C não vale na PH.
- WP-D0 (álgebra exata): na fase quebrada só N_ψ é conservado; zero cargas
  candidatas nas quatro variantes. A carga elétrica não é carga de Noether da
  eq. (1): o J·Σ completo a viola, e o salto 2↔3 do condensado a viola por ΔQ = 5/3.
- WP-D0′: com N_ψ total = 0 e matéria escura de antipsions,
  Ω_DM/Ω_b = 5,63·(m_DM/N_DM)/m_N, que exige m/N ≈ 0,89 GeV. Os 39 estados do
  censo dão −93% a −100%; a massa canônica do v85 (7,2×10¹⁰ GeV) erra para o outro
  lado. Saída consistente restante: excesso de antineutrinos de 1,88 por bárion,
  em conflito com μν/T = 0,123 de §13.5 se este for assimetria de número leptônico.

## Retificação de ontologia (19/10/2026)

As etapas WP-D0′ e WP-CP, e a afirmação "B+L = N_ψ/3 é conservado", supuseram
antimatéria feita de antipsions (ψ̄ de Dirac). O PH Foundations define a
antipartícula como inversão de índices k → (6−k) mod 6 (A.5) e atribui o excesso
a LLL ligar mais que RRR (J.1). Com essa ontologia, essas três conclusões caem;
WP-D0 (só tuplas de psions) continua válido.

## Resultados das etapas novas
- WP-J1: J·Σ é metade par e metade ímpar sob as duas conjugações (não é "ímpar").
  C_A acerta as cargas mas faz de e⁻ e p conjugados; C_B (Z2) não conjuga a carga.
- WP-M: com a eq. (10), ou com a eq. (12) lendo |m_J,R| como rótulo de excitação,
  (0,4,5) e (0,1,2) têm a mesma massa no mesmo nível. A função `_mass` de
  `P_yield_ratio.py` (|m_J,R| tirado da tupla) não segue o paper.
- WP-SIM: 12 corridas, controles OK; nenhuma produz hidrogênio, porque h e J·Σ são
  cegos ao nível e dão Δe = −Δp. Achado colateral: as cargas de V.1 somam +1.
  A regra de formação dependente do nível que faltaria é a mesma lacuna declarada
  em §8.2 (regra de níveis dos fragmentos, Teorema KS suspenso).

`simulador/` contém a página interativa (abre no navegador), o núcleo `core.js` e
`test.js`, que compara o JavaScript com o script Python (`node simulador/test.js`).

Pré-requisitos: Python 3, numpy, scipy; Node opcional para o teste do simulador.
