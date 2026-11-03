# WP-SPINF: a small, frozen family of spin rules for the nucleon magnetic form factors (pre-registration)

Why a family and not one rule: the minimal spin rule (WP-SPIN, rule R1) fails the neutron magnetic shape (WP-GM). Picking another rule after seeing that result is a fit. Here the whole family is declared BEFORE running, its size (7) is counted as a look-elsewhere factor, and a second stage on data not used for the selection (hyperon moments) is required for success.
Frozen inputs: extended-psion sizes r_0 = 0.413, r_1 = 0.442, r_4 = 0.691, r_5 = 0.340 fm, a = 0.940 fm (WP-FF3 B2); static hexagon; magnetisation model G_M^h(q) = kappa sum_i q_i <sigma_iz> f_i(q) j0(q a d_i) (WP-GM); charges Q9; shape reference G_D with 0.71 GeV^2; Q^2 = 0.1, 0.2, 0.3, 0.5, 0.8, 1.0 GeV^2.
State convention (all rules): H_s = sum_{i<j} J_ij s_i.s_j; the hadron is the equal-weight mixture over the ground manifold restricted to total S_z = +1/2 (a unique state when the ground doublet is non-degenerate). Rules with a spin-3/2 ground multiplet are excluded because the nucleon has spin 1/2.
The 7 rules (J_ij in units of J>0, antiferromagnetic unless stated):
 R1: all left-right pairs, equal weights (= WP-SPIN, known to fail the neutron shape).
 R2: all pairs, equal weights (Heisenberg triangle; mixture of the two doublets).
 R3: same-chirality pairs only, equal weights: that pair forms a singlet, the third psion carries the spin (scalar diquark).
 R4: same-chirality pairs only, ferromagnetic (J<0): the pair is a triplet; mixture over the degenerate S_z=+1/2 states.
 R5: all left-right pairs with weights 1/d_ij (chord length on the unit hexagon).
 R6: all pairs with weights 1/d_ij.
 R7: only the pair {2,3} (the bond of App. J.6.3), J>0; triples without both 2 and 3 have H_s = 0 (mixture over all S_z=+1/2 states).
Stage 1 (nucleons): a rule PASSES iff (i) |mu_n/mu_p - (-0.685)| <= 0.20 and (ii) the shape ratio G_M/G_M(0)/G_D is within 15% of 1 at all six points for BOTH p = (0,4,5) and n = (0,1,5).
Stage 2 (hyperons, not used in stage 1): for each rule passing stage 1: |mu_Lambda/mu_p - (-0.219)| <= 0.30 for Lambda = (1,2,3) and |mu_Sigma+/mu_p - 0.880| <= 0.30 for at least one of Sigma+ = (1,3,5), (2,3,4). (Wider tolerance because the hyperon moments are not shape-tested.)
Family success: at least one rule passes both stages. Family size 7: a passing rule is reported with 'look-elsewhere factor 7' and is NOT counted as confirmation unless it also passes the hyperon stage.
Observed values (PDG, quoted from memory; to be verified): mu_p = 2.793, mu_n = -1.913, mu_Lambda = -0.613, mu_Sigma+ = 2.458 (nuclear magnetons); ratios -0.685, -0.219, 0.880.
Anticipation (registered before running): no rule passes stage 1 (the neutron magnetisation weights q_i<sigma_i> either cancel and change sign, as in R1, or give the wrong sign of mu_n/mu_p); at most one rule passes the proton-only shape. Family success: NO.
Limits: static ring; no orbital part; the family is mine and small; Gaussian clouds; shape reference is the dipole (convention).
