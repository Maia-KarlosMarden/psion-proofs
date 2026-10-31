# WP-RULE: does a simple rule on the triple choose the geometry? (pre-registration)

Data (from WP-GEO, fixed): parsimony label = the simplest class compatible with the radius (order H < O < D;
H hexagon static, O octagon static, D dynamic): n (0,1,5):H; Sigma0 (0,2,4):H; Lambda (1,2,3):D; Sigma+a (1,3,5):O; Sigma+b (2,3,4):D.
The parsimony order is a CONVENTION; without it D is compatible with all five (radii alone need no geometry diversity).
Rule family: decision stumps on one binary feature of the triple, f=0 -> class a, f=1 -> class b (classes in {H,O,D}).
Features (14): contains psion j, j=0..5 (6); contains 2 or 3; contains both 2 and 3; sum k mod 6 == 0; #R psions >=2; #R == 1;
has an L-R pair; contains Z2-even psion count >=2 (k even); contains k mod 3 == 1 psion (always true by Z3 neutrality -> constant, kept as control).
Procedure: leave-one-out over the 5 composites: choose the stump (feature, a, b) with highest training accuracy (ties: first in the list order above,
classes ordered H,O,D), predict the held-out. Score = number correct of 5.
Baselines: (i) majority class of training; (ii) permutation test: all 120 permutations of the 5 labels, same procedure; p = fraction of permutations with score >= observed.
Success: observed score >= 4/5 AND permutation p <= 0.05. Look-elsewhere: only the family above was tried (frozen before running).
Anticipation (registered before running): FAILS. Score <= 3/5 and p > 0.2: five points and three classes do not constrain a stump.
