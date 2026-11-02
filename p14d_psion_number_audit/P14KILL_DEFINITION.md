# WP-KILL: kill criteria for "the subatomic zoology is made of psion composites" and an attempt to falsify them (pre-registration)

CLAIM (as the paper states it): particles are Z6 composites of psions; the unconfined triples carry the charges of Q9=(1,-5,-5,10,4,4)/9;
the ground state of a hadron-like triple is a configuration of its psions on a ring (static hexagon at the natural sites); masses sit on a cascade ladder m(l)=Lambda_meta exp(-l/b_e), ratio 1/34.6 per rung,
anchored by the electron at l=10 (m=0.5089 MeV); the nucleon scale is identified, at leading order, with l=8 (610 MeV).
Auxiliary models tried in this audit (dynamic rings, hybrid, geometry per composite, minimal spin, K.10 vacua) are NOT part of the claim; their failure is recorded elsewhere and does not count here.

KILL CRITERIA (parameter-free; no quantity is fitted):
 K1 (radius of the Sigma-, static hexagon, no parameters). Prediction R_pred = <r2>_Sigma-/<r2>_p from the ring at natural sites (same code convention as WP-GEO).
    Data (SELEX 2001, quoted from a secondary source): <r2>_Sigma- = 0.61 +- 0.21 fm^2; proton <r2>=0.8414^2=0.7079 fm^2 -> R_obs = 0.862 +- 0.297. Control: neutron R_pred vs R_obs=-0.1161/0.7079=-0.164 +- 0.003 (+- 0.0164 as used in the audit window).
    KILL if |R_pred - R_obs| > 2 sigma (0.594) for the Sigma-. SUPPORT if <= 1 sigma (0.297).
 K2 (cascade ladder alignment). l_eff(m) = 10 - ln(m/0.5089 MeV)/ln(34.6) with ln(34.6)=3.5446 (paper values). A composite on a rung has l_eff integer.
    Sample: 25 light and heavy hadron ground states (PDG masses quoted from memory, to be verified): pi+-, pi0, K+-, K0, eta, rho, omega, eta', phi, p, n, Lambda, Sigma+, Sigma0, Sigma-, Xi0, Xi-, Delta(1232), Sigma*(1385), Xi*(1530), Omega-, D0, J/psi, B+-, Upsilon(1S).
    Statistic: N_in = number with |l_eff - round(l_eff)| <= 0.15 (0.15 chosen because the paper's own nucleon identification is 0.12 rung off). Null: grid-offset permutation (shift all l_eff by a uniform offset in [0,1), 10^5 draws); p = P(N_null >= N_in).
    KILL if p > 0.20; SUPPORT if p <= 0.05; otherwise inconclusive.
 K3 (neutron-proton mass sign, structural). In the ring the Coulomb sum B gives E_n - E_p = -0.43 kappa (computed earlier). With kappa>0 (physical repulsion) m_n < m_p, observed m_n - m_p = +1.293 MeV. Recorded as FAIL under the stated premise (not rerun); counted as a failed criterion of the claim only if the paper attributes the n-p splitting to the ring Coulomb energy: it does not state so. Therefore K3 is reported but NOT counted.
 K4 (minimal spin, magnetic moments): reported from WP-SPIN, NOT counted (auxiliary postulate).
 K5 (counting of states): NOT decidable, because the cascade level l is an independent quantum number (the same triple appears at several levels) so 8 triples do not bound the number of states. Reported as 'no test'.
VERDICT RULE (frozen): the claim is FALSIFIED in the registered form if K1 and K2 both come out KILL; SUPPORTED if both SUPPORT; otherwise UNDETERMINED (one of the two decisive tests failed or is inconclusive: the claim is then 'not refuted, not supported').
Anticipation (registered before running): K1 within 2 sigma (wide SELEX error) but not within 1 sigma; K2 KILL (p > 0.2): hadron masses cluster near 0.5-1 GeV because of the QCD scale, not because of the electron-anchored grid. Overall UNDETERMINED.
Look-elsewhere: only K1 and K2 are decisive; both frozen here.
