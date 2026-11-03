# WP-JNT: joint fit of electric and magnetic nucleon form factors with extended psions (pre-registration)

Question: can ONE set of psion sizes (r_0, r_1, r_4, r_5) accommodate the proton and neutron, electric and magnetic, at once, for a fixed spin rule? 
Rules fixed beforehand (the two that passed the proton magnetic shape in WP-SPINF): R1 (all left-right pairs, equal weights; = WP-SPIN) and R4 (like-chirality pairs only, ferromagnetic). No other rule is tried.
Model: static hexagon, Gaussian psion clouds, ring radius a fixed by the proton slope (0.8414 fm); electric terms exactly as WP-FF3; magnetic G_M^h(q) = sum_i q_i <sigma_iz> f_i(q) j0(q a d_i) with the spin expectations of the chosen rule (state convention of WP-SPINF); shape G_M/G_M(0)/G_D at Q^2 = 0.1, 0.2, 0.3, 0.5, 0.8, 1.0 GeV^2.
Loss L = max of six terms: (1) electric proton shape: max|F_p/G_D - 1|/0.15; (2) electric neutron shape: max|F_n - G_E^n|/(0.15 max|G_E^n|); (3) R_n: |R_n + 0.164|/0.0164; (4) magnetic proton shape: max|G_M^p/G_M^p(0)/G_D - 1|/0.15; (5) magnetic neutron shape: max|G_M^n/G_M^n(0)/G_D - 1|/0.15; (6) moment ratio: |mu_n/mu_p + 0.685|/0.20 (fixed by the rule, independent of the sizes). Viable iff L < 1.
Search: differential evolution, popsize 40, maxiter 300, ten seeds (91..100) per rule, bounds r_i in [0, 0.9] fm, polish on.
Decision: SUPPORT-OF-CONSISTENCY if some rule has L < 1 (descriptive: 4 parameters for ~24 numbers; the sizes are then reported; nothing is predicted); otherwise the rule is declared incompatible with the extended-psion hexagon.
Anticipation (registered before running): NOT viable for either rule; best L between 1.5 and 3.0; the binding term is the magnetic neutron shape.
Limits: static ring, no orbital magnetisation, no recoil; dipole shape reference (a convention) for both G_M; tolerances mine; PDG ratio from memory.
