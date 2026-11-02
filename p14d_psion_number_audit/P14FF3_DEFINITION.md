# WP-FF3: static hexagon with extended psions (pre-registration)

Motivation (user): psions have volume; WP-FF and WP-FF2 used point psions. The paper (v86) states NO psion size; r_psion is therefore a new quantity, not derived.
Model: static hexagon (p=(0,4,5), n=(0,1,5), as in L.4.7), each psion i with a Gaussian charge cloud of rms radius r_i (3D): form factor f_i(q)=exp(-q^2 r_i^2/6). Hadron:
 F_h(q) = sum_i q_i f_i(q) j0(q a d_i),  with d_i the unit-ring distance of psion i from the unweighted centroid of the three.
 Slope: <r2>_h = [ sum_i q_i (a^2 d_i^2 + r_i^2) ] / Q_h for Q_h != 0 (p), and the unnormalised sum for the neutron (as before).
 The ring radius a is fixed by the proton slope: <r2>_p = 0.8414^2 = 0.7079 fm^2.
Version A (common size): r_i = r_psi for all psions; scan r_psi in [0, 0.9] fm (91 points). Neutron slope R_n = <r2>_n/<r2>_p must lie in -0.164 +- 0.0164 (z-score); shape comparison as in WP-FF (dipole G_D and Galster G_E^n at Q^2 = 0.1,0.2,0.3,0.5,0.8,1.0 GeV^2).
Version B (per-psion sizes): r_0, r_1, r_4, r_5 free in [0, 0.9] fm (the four labels used by p and n), five DE starts (seeds 72..76, popsize 15, maxiter 40); same loss. This is a fit with 4 parameters to 12 numbers, not a prediction.
Loss L = max( max_Q |F_p/G_D - 1|/0.15 , max_Q |F_n - G_E^n|/(0.15 max|G_E^n|) , |R_n + 0.164|/0.0164 ). Viable iff L < 1.
Decision: version A: SUPPORT if some r_psi has L < 1, otherwise KILL. Version B: reported as DESCRIPTIVE; if viable, report the sizes and whether they are of one order of magnitude (max/min <= 3, excluding any r_i < 0.02 fm).
Anticipation (registered before running): A: KILL; the neutron slope limits r_psi to <= 0.17 fm (R_n=-0.154(1 - r_psi^2/0.708) must stay <= -0.1476) and such a size changes the form factor by < 6% at Q^2 = 0.5 GeV^2, so the 38% deficit remains. B: viable L < 1 is possible by choosing unequal sizes; I anticipate that it is reached only with sizes differing by more than a factor 3, i.e. a fine-tuned description.
Limits: Gaussian clouds are a choice; no motion; no recoil; dipole/Galster parametrisations (WP-FF caveats).
