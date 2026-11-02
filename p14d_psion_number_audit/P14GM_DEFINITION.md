# WP-GM: magnetic form factors with the frozen extended-psion sizes (pre-registration)

Out-of-sample use of WP-FF3 (B2 optimum, frozen): r_0 = 0.413, r_1 = 0.442, r_4 = 0.691, r_5 = 0.340 fm, ring radius a = 0.940 fm; psion Gaussian clouds f_i(q) = exp(-q^2 r_i^2/6). These were fitted to the ELECTRIC shape and R_n only; no magnetic quantity was used.
Spin (WP-SPIN, frozen): H_s = J sum_{L-R pairs} s_i.s_j with J>0 (ground multiplet J=1/2, S_z=+1/2); moment operator mu = kappa sum_i q_i sigma_iz (kappa cancels); spatial and spin parts factorise (static ring).
Model magnetic form factor: G_M^h(q) = kappa sum_i q_i <sigma_iz> f_i(q) j0(q a d_i) (spin magnetisation sitting on the psion clouds; no orbital or tensor term). Moments: G_M^h(0)=kappa sum_i q_i <sigma_iz>; the ratio mu_n/mu_p is as in WP-SPIN (-0.800, observed -0.685).
Shape comparison (same points as WP-FF): ratio G_M^h(Q^2)/G_M^h(0) against G_D(Q^2)=(1+Q^2/0.71 GeV^2)^-2, for p and n, Q^2 = 0.1, 0.2, 0.3, 0.5, 0.8, 1.0 GeV^2. The measured G_M^p/(mu_p G_D) and G_M^n/(mu_n G_D) follow the dipole only approximately (arXiv:1302.4604 says neither is well represented by a dipole from t=0); tolerance 15% as in WP-FF, with that caveat.
Decision: SUPPORT if |ratio/G_D - 1| <= 0.15 at all six points for both p and n; KILL if p or n fails at >= 3 of the six points; otherwise INCONCLUSIVE.
Anticipation (registered before running): KILL. The magnetisation weights q_i<sigma_iz> differ from the charge weights q_i that fixed the sizes (for the neutron they can nearly cancel), so the magnetic shape is not given by the electric fit.
Limits: static ring, no orbital magnetisation, no relativistic or recoil corrections; the spin postulate failed on the hyperon moments (WP-SPIN), so this is the same spin model on the nucleon only.
