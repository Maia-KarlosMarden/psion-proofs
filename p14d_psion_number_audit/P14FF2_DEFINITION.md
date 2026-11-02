# WP-FF2: form factors with a delocalised (dynamic 12-site) ring (pre-registration)

Follows WP-FF (static hexagon fails the shape of the elastic form factors). Other version tested here: p and n as the ground state of the DYN12 Hamiltonian H=diag(A+kB+lC)-tT-sS (same code as WP-GEO/ESC; 1320 placements on 12 sites), with the form factor
F_h(q) = sum_placements |psi(placement)|^2 sum_i q_i j0(q |x_i - c|), c the unweighted centroid of the three psions in that placement; orientation average as before.
Parameters (t,s,k,l) shared by p and n; bounds t,s in [0,6], k,l in [0,20]. The ring radius a is fixed per parameter point by the proton slope: a = 0.8414 fm / sqrt(<R2>_p) (parameter-free given (t,s,k,l)).
Data and comparison points: exactly as WP-FF (dipole G_D with 0.71 GeV^2 and Galster G_E^n; Q^2 = 0.1, 0.2, 0.3, 0.5, 0.8, 1.0 GeV^2; hbar c = 0.19733 GeV fm); additionally the neutron slope R_n=<r2>_n/<r2>_p must lie in -0.164 +- 0.0164 (z-score).
Loss L = max( max_Q |F_p/G_D - 1| / 0.15 , max_Q |F_n - G_E^n| / (0.15 max|G_E^n|) , |R_n + 0.164| / 0.0164 ). Viable iff L < 1.
Search: differential evolution (popsize 15, maxiter 40, seeds 72..76; five independent starts); ARPACK failures count as invalid points (loss 50).
Reference: the static hexagon has L = 2.55 (proton at Q^2=0.5) vs 0.38/0.15, and neutron 0.78/0.15 = 5.2, so L_static = 5.2.
Decision: SUPPORT if >=2 of 5 starts are viable and the viable optima also pass the same shape test at two EXTRA points (Q^2 = 0.4 and 0.6 GeV^2, not used in the loss) within the same tolerances; INCONCLUSIVE if the viable optima fail the extra points; KILL if no start is viable.
Anticipation (registered before running): KILL (no start with L<1); best L between 1.5 and 4. Delocalisation smooths the structure but cannot give the proton shape and the neutron shape and R_n at once with 4 parameters.
Limits: this is a fit with four free parameters to twelve numbers, not a prediction; viability would be descriptive; the hybrid HYB3 (p, n static) is NOT this version: with p, n static it reduces to WP-FF and fails.
