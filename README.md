# Psion Hypothesis — v8.4

**A Framework for the Origin of Matter, Antimatter, Dark Matter, and Gravity**

Karlos Marden Maia · Independent Researcher, Marietta, Georgia
`karlosmardenmaia@gmail.com`

---

## Overview

The Psion Hypothesis (PH) derives Standard Model constants and cosmological
parameters from a single fermionic field ψ with Z₆ = Z₃ × Z₂ internal
symmetry and spin j = 5/2.

**The framework has two empirical (SC) inputs:** Newton's constant **G_N**
(the cascade anchor) and **Δm²₂₁** (the neutrino-sector scale, §11). Each has
its own, independent closure programme targeting zero remaining inputs —
neither is closed yet:

```
m_e = Λ_meta · exp(−N_e/b_e)          [Λ_meta fixed from G_N via Sakharov
                                        induced gravity; N_e=10, b_e=325/1152
                                        are exact Z6 theorems, T1]
```

- **P2/P12** ([T2]): solve the fixed-point condition Φ(Λ)=Λ of the two-loop
  NJL kernel → closes G_N, SC inputs 2→1.
- **P26-C** ([T4], not promoted): candidate closed form for Δm²₂₁ → would
  close SC inputs 1→0. Registered but not adopted.

### Key Results (v8.4, PH_Foundations §11–12)

| Observable | Predicted | Observed | Error | Status |
|---|---|---|---|---|
| G_N (GeV⁻²) | 6.70883×10⁻³⁹ | 6.70883×10⁻³⁹ | 0.0% | [SC] input |
| G_N (GeV⁻²), cascade-recovered | 6.653×10⁻³⁹ | 6.709×10⁻³⁹ | 0.83% | [T2] |
| m_e (MeV) | 0.5089 | 0.511 | 0.42% | cascade output |
| Δm²₂₁ (eV²) | 7.530×10⁻⁵ | 7.530×10⁻⁵ | 0.0% | [SC] input |
| ρ_Λ (GeV⁴) | 2.492×10⁻⁴⁷ | 2.515×10⁻⁴⁷ | 0.93% | [T2+] |
| ρ_Λ·G_N² | 1.121×10⁻¹²³ | 1.132×10⁻¹²³ | 0.93% | [T2+] |
| η_B | 6.02×10⁻¹⁰ | 6.00×10⁻¹⁰ (6.13×10⁻¹⁰ base-chain) | 0.4% / 1.7% | [T2] |
| Ω_DM/Ω_b | 5.43 (direct) / 5.64 (inference) | 5.364 | 1.2% / 5.1% | [T2] |
| N_eff | 3.054 | 3.04 ± 0.18 | <1σ | [T2] |
| N_ν | 3 | 3 | exact | [T1] |
| N_g | 8 | 8 | exact | [T1] |
| N_gen (charged leptons) | 3 (kinematic bound, JR=7/2>jmax) | 3 | — | [T1] |
| DM stable states | 39 (N∈{9,12,15,18,21,27}, term. to N=48) | — | — | [T1] |

**17 falsifiable predictions (F1–F17)** are registered in the paper's Table
4, spanning four timescales, with T1 hard falsifiers on N_ν, N_g, and N_gen.
None of F1–F17 depends on the still-open P8b lepton-mass refinement (see
Appendix L.4 of the paper).

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
│   ├── PH_Foundations_v84.tex        ← LaTeX source (canonical; see note below)
│   └── PH_Foundations_v84.pdf        ← compiled PDF, 43 pages
├── proofs/
│   ├── scripts/
│   │   ├── cascade_anchor.py         ← SC input declaration
│   │   └── ...
│   ├── tests/
│   │   └── ...
│   ├── proofs/
│   │   ├── P_dark_matter.py          ← 39-state DM census (2 independent algorithms)
│   │   ├── P_dark_energy.py
│   │   ├── P_eps_psi.py              ← chirality asymmetry ε_ψ = 25/169
│   │   ├── P_counting_theorems.py    ← N_ν=3, N_g=8, N_gen=3
│   │   ├── P_charge_algebra.py
│   │   ├── P_body_number_law.py
│   │   ├── P_wall_tensions.py
│   │   ├── P_yield_ratio.py
│   │   ├── P_kappa_reduction.py
│   │   └── P_predictions.py
│   └── docs/
│       ├── CHARGE_CONVENTIONS.md     ← proton/neutron/exotic charge audit;
│       │                                conjugation impossibility proof (§3b);
│       │                                adopted resolution (§3c)
│       ├── SECTION_12_NOTE.md        ← companion note on the C/CP reading
│       │                                used in the paper's Baryogenesis §13
│       ├── PSION_DIAGRAM_TEXTBOOK.md ← worked PD-protocol examples; flags an
│       │                                open e+/proton tuple-collision item
│       ├── QFT_PIPELINE.md           ← baryogenesis yield-ratio pipeline
│       ├── PROGRESS_AND_ROADMAP.md   ← version history, closure roadmap
│       └── DARK_MATTER_PROOFS.md
└── docs/
    └── CHANGES_v83.md
```

**Note on versioning:** the canonical source is the `.tex` file. An earlier
`.docx` export of this paper diverged from the `.tex` in the Baryogenesis
section wording; the `.tex` is authoritative. If a `.docx` is regenerated for
distribution, regenerate it *from* the current `.tex`, not edited
independently — see `PROGRESS_AND_ROADMAP.md` for the record of this issue.

---

## v8.4 Audit Trail (this revision)

v8.4 closed a self-directed audit covering charge conventions, the
falsifiability registry, and the SC-input count. Full detail is in
`CHARGE_CONVENTIONS.md` and the paper's Appendix L.4–L.5. Summary:

- **Composite identification corrected:** proton = LRR (0,4,5), Q=+1;
  neutron = LLR (0,1,5), Q=0; RRR (3,4,5) is a separate Q=+2 exotic state.
  Both proton and neutron charges hold exactly under Theorem V.1 alone.
- **State-wise charge conjugation does not hold under Theorem V.1** (T1,
  independent of particle labelling — a property of the six Q_k values
  themselves). Two candidate repairs were closed by exhaustive search;
  conjugation is retained only as an aggregate identity over the composite
  spectrum. This does not affect η_B, which never invokes state-wise
  Q_k-conjugation (paper §13.1, footnote).
- **μ⁻/τ⁻ are not N=3 tuples** — Theorem 6.2 (unique Q=−1 three-body state)
  forecloses that reading. Their mechanism is the Leptonic Transfer Matrix
  (P8b, paper App. A.2.5); the formal integration step (Passo A.3) remains
  [T2], not yet [T1]. **None of the 17 falsifiable predictions depend on
  this closure.**
- **SC-input count corrected from 1 to 2** (G_N and Δm²₂₁ are independent
  inputs with independent, separately-tracked closure candidates).
- **Falsifiable-prediction count corrected to 17** (a stale "16" in one
  table caption did not match the 17 rows actually listed, nor the abstract).
- **Open item, not resolved:** `PSION_DIAGRAM_TEXTBOOK.md` derives e⁺ =
  (0,4,5) by direct Z6 conjugation of the electron tuple — the same tuple
  Appendix H assigns to the proton. Both are charge-consistent (Q=+1) but
  nothing in the framework currently distinguishes a lepton antiparticle
  from a baryon at the same tuple. Flagged, not adjudicated.

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
| P8b (Passo A.3) | Formal integration of Σ₂ = b_e×15/8 | [T2] | lepton-mass errors 5–9% → <1% (no F-prediction depends on this) |
| P23 | Baryogenesis closure — two components still open | [T2+] | η_B, ε_ψ → T1 |
| P1 | b_e = 325/1152 from 2PI kernel | [T2] | mass ratios Class B |
| P-4π/3 | Prefactor 4π/3 in Ω_DM/Ω_b | [T2+] | Ω_DM → T1 |
| — | e⁺/proton tuple collision at (0,4,5) | open, unresolved | see Audit Trail above |

---

## Hard Falsifiers (T1)

The following experimental results would **immediately and irrevocably**
falsify the PH without further calculation:

- Detection of a **4th charged lepton**
- Detection of a **4th active neutrino**
- Detection of a **9th gluon**
- Any stable neutral Z₆ composite at N=4 or N>27 (outside {9,12,15,18,21,27})

---

## Citation

```bibtex
@misc{maia2026ph,
  author    = {Maia, Karlos Marden},
  title     = {The Psion Hypothesis: A Framework for the Origin of
               Matter, Antimatter, Dark Matter, and Gravity},
  year      = {2026},
  note      = {PH\_Foundations v84, \url{https://github.com/karlosmardenmaia/psion-proofs}}
}
```

---

## License

© 2026 Karlos Marden Maia. All rights reserved.
The companion verification repository `psion-proofs` is open-source (MIT).
