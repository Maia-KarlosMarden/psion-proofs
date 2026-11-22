# Psion Hypothesis — v8.6

**A Framework for the Origin of Matter, Antimatter, Dark Matter, and Gravity**

Karlos Marden Maia

---

## Overview

The Psion Hypothesis (PH) derives the Standard Model's structural constants
and the mass–gravity–cosmology observables (CKM/PMNS mixing and individual
quark masses are not derived) from a single fermionic field ψ with
Z₆ = Z₃ × Z₂ internal symmetry and spin j = 5/2.

**The framework has two empirical (SC) inputs:** Newton's constant **G_N**
(the cascade anchor) and **Δm²₂₁** (the neutrino-sector scale, §11). Each has
its own, independent closure programme targeting zero remaining inputs —
neither is closed yet:

```
m_e = Λ_meta · exp(−N_e/b_e)          [Λ_meta fixed from G_N via Sakharov
                                        induced gravity; b_e = 325/1152 is
                                        fixed by Z6 (T1); N_e = 1+8+1 = 10 is
                                        exact given N_g = 8 (P-Ng)]
```

- **P2/P12** ([T2]): solve the fixed-point condition Φ(Λ)=Λ of the two-loop
  NJL kernel → closes G_N, SC inputs 2→1.
- **P26-C** ([T4], not promoted): candidate closed form for Δm²₂₁ → would
  close SC inputs 1→0. Registered but not adopted.

### Key Results (v8.6; PH_Foundations §11–14)

With G_N as the anchor, the mass–gravity–cosmology sector has **two
independent numerical residuals**: m_e (0.42%) and the invariant ρ_Λ G_N²
(0.93%). ρ_Λ and G_N m_e² are re-expressions of these.

| Observable | Predicted | Observed | Error | Status |
|---|---|---|---|---|
| G_N (GeV⁻²) | 6.70883×10⁻³⁹ | 6.70883×10⁻³⁹ | 0.0% | [SC] input |
| m_e (MeV) | 0.5089 | 0.511 | 0.42% | cascade output |
| Δm²₂₁ (eV²) | 7.530×10⁻⁵ | 7.530×10⁻⁵ | 0.0% | [SC] input |
| ρ_Λ (GeV⁴) | 2.492×10⁻⁴⁷ | 2.515×10⁻⁴⁷ | 0.93% | [T2+] |
| ρ_Λ·G_N² | 1.121×10⁻¹²³ | 1.132×10⁻¹²³ | 0.93% | [T2+], G_N-independent |
| η_B | 6.02×10⁻¹⁰ | 6.13×10⁻¹⁰ (base-chain) | 1.7% | consistency relation, given ξ* (P14-ξ*) |
| Ω_DM/Ω_b | 5.43 (direct) / 5.64 (inference) | 5.364 | 1.2% / 5.1% | direct route tied to η_B via ξ* |
| N_eff | 3.064 | 3.04 ± 0.18 | <1σ | [T2]; μ_ν/T in tension with BBN (T-BBN) |
| N_ν | 3 | 3 | exact | [T1] |
| N_g | 8 | 8 | — | [T3], conditional on SU(3) emergence (P-Ng) |
| N_gen (charged leptons) | 3 (kinematic bound, JR=7/2>jmax) | 3 | — | [T1] |
| DM stable states | 39 (N∈{9,12,15,18,21,27}, term. to N=48) | — | — | [T1] |
| m_DM | 2Λ_meta/35 = 7.2×10¹⁰ GeV | — | — | [T2] |

The dark matter of the PH is cold and collisionless at all astrophysical
scales (§14.5). The RGD kernel has been falsified by mock controls.

**18 falsifiable predictions (F1–F17, with F5 split into F5a/F5b)** are
registered in the paper's Table 4, spanning four timescales. None depends on
the still-open P8b lepton-mass refinement (Appendix L.4 of the paper).

### Epistemic Classification

| Label | Meaning |
|---|---|
| **T1** | Exact algebraic theorem — no numerical verification required |
| **T2** | Structurally derived and numerically verified |
| **T2+** | One explicit algebraic step from T1 |
| **T3** | Pipeline defined, calculation open |
| **T4** | Conjectural: mechanism identified but not derived |

---

## Repository Structure

```
ph_v83_repo/
├── paper/
│   ├── PH_Foundations_v86.tex        ← LaTeX source, current version (canonical)
│   ├── PH_Foundations_v86.pdf        ← compiled PDF, 49 pages
│   ├── PH_Foundations_v85.tex        ← first published version
│   ├── PH_Foundations_v85.pdf
│   ├── PH_Foundations_v84.tex        ← pre-publication draft
│   └── PH_Foundations_v84.pdf
├── proofs/
│   ├── scripts/
│   │   ├── cascade_anchor.py         ← SC input declaration
│   │   └── ...
│   ├── tests/
│   │   └── ...
│   ├── proofs/
│   │   ├── P_dark_matter.py          ← 39-state DM census (fragmentation criterion, 2 algorithms)
│   │   ├── P_kinematic_stability.py  ← inequality (11) on level partitions; KS scope guard
│   │   ├── P_pmns_structure.py       ← PH flavour basis (PMNS programme, Stage 0)
│   │   ├── P_dark_energy.py
│   │   ├── P_eps_psi.py              ← chirality asymmetry ε_ψ = 25/169
│   │   ├── P_counting_theorems.py    ← N_ν=3, N_gen=3; N_g=8 given SU(3) (P-Ng)
│   │   ├── P_charge_algebra.py
│   │   ├── P_body_number_law.py      ← RETRACTED for the DM mass (kept for the record)
│   │   ├── P_wall_tensions.py
│   │   ├── P_yield_ratio.py
│   │   ├── P_kappa_reduction.py
│   │   └── P_predictions.py
│   └── docs/
│       ├── CHARGE_CONVENTIONS.md     ← proton/neutron/exotic charge audit
│       ├── SECTION_12_NOTE.md        ← companion note on the C/CP reading (§13)
│       ├── PSION_DIAGRAM_TEXTBOOK.md ← worked PD-protocol examples
│       ├── QFT_PIPELINE.md           ← baryogenesis yield-ratio pipeline
│       ├── PROGRESS_AND_ROADMAP.md   ← version history, closure roadmap
│       ├── DARK_MATTER_PROOFS.md
│       └── PMNS_PROGRAMME.md         ← pre-registered PMNS model log; F3 status
└── docs/
    ├── CHANGES_v83.md
    ├── CHANGES_v84.md
    ├── CHANGES_v85.md
    └── CHANGES_v86.md
```

**Note on versioning:** the canonical source is the `.tex` file. If a `.docx`
is generated for distribution, generate it *from* the current `.tex`.

---

## v8.6 (current version)

Second published version. An external audit recomputed every verifiable number
from the declared inputs; the core arithmetic reproduces to ≤ 0.1% and no core
number changed. Full record in `docs/CHANGES_v86.md`. Summary:

- **Success ledger:** the "G_N recovered at 0.83%" row is twice the m_e
  residual, not an independent test; two independent residuals remain.
- **Baryogenesis:** η_B and the direct Ω_DM/Ω_b route form one consistency
  relation until ξ* is derived independently (P14-ξ*).
- **Gluon count:** N_g = 8 (and N_e = 10, F9) conditional on the emergence of
  SU(3) colour (P-Ng); SU(6) shares the centre Z6.
- **Neutrinos:** N_eff = 3.064 (standard coefficient); BBN tension T-BBN
  registered; Δm²₃₁ comparator corrected.
- **Consistency:** CEC vacuum energy vs ρ_Λ registered (PC3); ANITA
  reaction removed pending a charge-conserving derivation; Table 12 rebuilt.
- **Dark matter in halos:** new §14.5 (cold, collisionless); RGD falsified.
- **Falsifiable predictions:** counted as 18.

## Earlier versions

- **v8.5** (first published version): DM.1 census proof fixed; canonical DM
  mass at the Landau pole; Theorem KS physical application suspended; F5 split
  into F5a/F5b with an ex-ante kill criterion; F3 conditional on P-PMNS. See
  `docs/CHANGES_v85.md`.
- **v8.4**: charge-convention audit (proton LRR (0,4,5), neutron LLR (0,1,5));
  state-wise charge conjugation retained only as an aggregate identity;
  μ⁻/τ⁻ via the Leptonic Transfer Matrix (P8b); SC-input count 2. See
  `docs/CHANGES_v84.md` and `proofs/docs/CHARGE_CONVENTIONS.md`.

---

## Running the Proofs

```bash
pip install -r proofs/requirements.txt
cd proofs

# Full suite
bash scripts/run_all_proofs.sh

# Single SC input declaration
python scripts/cascade_anchor.py

# 39-state dark matter census (two independent algorithms)
python -m proofs.P_dark_matter
```

---

## Open Problems

| ID | Problem | Status | Impact |
|---|---|---|---|
| P2/P12 | Derive Λ_meta from Φ(Λ)=Λ without an external G_N input | [T2] | SC inputs 2→1 |
| P26-C | Closed-form candidate for Δm²₂₁ (registered, not promoted) | [T4] | SC inputs 1→0 |
| P-Ng | Derive SU(3) colour from the condensate | open | N_g = 8, N_e = 10, F9 |
| P14-ξ* | Independent derivation of the supercooling parameter ξ* | [T3] | η_B and Ω_DM/Ω_b become predictions |
| P-ε_conf | Re-derive the ε_conf = 4/9 count | open | η_B chain |
| PC3 | CEC vacuum energy vs ρ_Λ | open | Prop. J.3, [G2] |
| T-BBN | μ_ν/T = 0.123 vs BBN bounds; write out the δL → μ_ν/T map | tension | F12 |
| P8b (Passo A.3) | Formal integration of Σ₂ = b_e×15/8 | [T2] | lepton-mass errors 5–9% → <1% (no F-prediction depends on this) |
| P1 | b_e = 325/1152 from 2PI kernel | [T2] | mass ratios Class B |
| P-4π/3 | Prefactor 4π/3 in Ω_DM/Ω_b | [T2+] | Ω_DM → T1 |
| — | e⁺/proton tuple collision at (0,4,5) | open | see `PSION_DIAGRAM_TEXTBOOK.md` |

---

## Hard Falsifiers (T1)

The following experimental results would falsify the PH without further
calculation:

- Detection of a **4th charged lepton**
- Detection of a **4th active neutrino**
- Any stable neutral Z₆ composite at N=4 or N>27 (outside {9,12,15,18,21,27})

A **9th gluon** falsifies the PH conditionally on P-Ng (SU(3) colour).

---

## Citation

```bibtex
@misc{maia2026ph,
  author    = {Maia, Karlos Marden},
  title     = {The Psion Hypothesis: A Framework for the Origin of
               Matter, Antimatter, Dark Matter, and Gravity},
  year      = {2026},
  note      = {PH\_Foundations v86, \url{https://github.com/Maia-KarlosMarden/psion-proofs}}
}
```

---

## License

© 2026 Karlos Marden Maia. All rights reserved.
The companion verification repository `psion-proofs` is open-source (MIT).
