# WP-MIX4: corrected rerun (static proton) + out-of-sample for mixed geometries (pre-registration)
Correction: radii ratios R_h = <r2>_h / <r2>_p with the STATIC proton (home placement, t=s=0). Everything else as WP-MIX stage 2
(loss, DIRS, DE popsize 15 maxiter 40, bounds). Models MIXa, MIXb, HYB. Starts: DE seeds 72..81 (10 per model); an optimum counts if loss<1.
Geometry of the out-of-sample composites is frozen BEFORE running with the HYB rule (dynamic iff contains psion 2 or 3): Sigma- (0,1,2) and X=(3,4,5) are D12 dynamic
at the optimum parameters (the static alternative is not used).
Out-of-sample criteria (unchanged from WP-OOS): P-A: 0.65<=R_Sm<=1.07 and m(Sigma-)-m(n) in [181,335] MeV with scale u=176.1/D (D=M_L-M_n) for at least one V0-approved direction;
P-B: m(X)-m(n) in [205,381] MeV for the same direction (both must hold for the same direction). V0 = radii + masses as in the model.
Score per model: fraction of viable optima (loss<1) passing P-A&P-B for at least one direction. Success for a geometry: viable optima >=3 and fraction >=0.6 (the
HYB reference in WP-OOS was 1 of 3).
Anticipation (registered before running): the corrected rerun keeps all three models viable (static p changes loss by <0.15), but no model reaches fraction 0.6 (HYB ~0.3, MIX <=0.3).
