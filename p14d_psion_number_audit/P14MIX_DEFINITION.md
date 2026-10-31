# WP-MIX: mass as second constraint on per-composite geometry (pre-registration; descriptive)

Question: if the composites take the geometries suggested by the radii (WP-GEO), does the mass criterion V0 still hold?
Mass functional as before: M_h = a<A>+b<B>+g<C> with the 4000 natural +/- directions of HYB3 (DIRS); V0 = D=M_L-M_n>0,
|(M_S0-M_L)/D/0.437-1|<0.3, |(M_Sp-M_p)/D/1.425-1|<0.3, q=D/max|dir|>=0.336. Same A,B,C normalisation (unit ring radius) in every geometry.
Models (p, n always hexagon static; <A,B,C> at the home placement, from the 12-site code at t=s=0,k=l=0):
 MIXa: Sigma0 hexagon static; Lambda D12 dynamic; Sigma+ = (1,3,5) octagon static (8-site code at t=s=0).
 MIXb: Sigma0 hexagon static; Lambda and Sigma+ = (2,3,4) D12 dynamic with shared (t,s,k,l).
 HYB (reference, rule of HYB3): Sigma0, Lambda and either Sigma+ candidate D12 dynamic, p,n static.
Dynamic parameters: N=4000 LHS points (seed 71; t,s log-uniform [0.01,6] with 5% zeros; k,l uniform [0,20]).
Radius criteria (strict): Lambda R in [0.05,0.25]; Sigma+ R in [0.80,1.15]; Sigma0 (HYB only) |R|<=0.15.
Output per model: points with >=1 direction satisfying radii AND V0; pairs (point,direction).
Success for the per-composite-geometry idea: MIXa or MIXb has >=1 viable point. Reference: HYB count.
Anticipation (registered before running): MIXa and MIXb have ZERO viable points while HYB has some. Reason: with Sigma0 static the ratio
(M_S0-M_L)/D >0 requires Sigma0 heavier than the dynamic Lambda, but delocalisation raises <A>, so mass forces Sigma0 to be dynamic too.
Look-elsewhere: no extra freedom added beyond the three models above.
