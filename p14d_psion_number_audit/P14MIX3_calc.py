import numpy as np
import P14MIX2_calc as W
OPT={"HYB":[0.462,0.497,0.021,0.018],"MIXa":[1.064,0.198,0.0,0.139],"MIXb":[0.442,0.389,0.025,0.04]}
rng=np.random.default_rng(73)
for m,x in OPT.items():
    x=np.array(x); L=[]
    for i in range(50):
        y=x*(1+0.1*rng.uniform(-1,1,4)); L.append(W.loss(y,m))
    L=np.array(L); print(f"{m}: loss nominal {W.loss(x,m):.3f}; perturbados: fracao <1 = {(L<1).mean():.2f}; mediana {np.median(L):.3f}; max {L.max():.3f}")
