import numpy as np
from scipy.optimize import differential_evolution
import P14FF3_calc as F
f=lambda x: F.Lof({0:x[0],1:x[1],4:x[2],5:x[3]})
res=[]
for sd in range(81,91):
    r=differential_evolution(f,[(0,0.9)]*4,popsize=40,maxiter=300,seed=sd,tol=1e-10,updating='immediate',polish=True)
    o=F.evalr({0:r.x[0],1:r.x[1],4:r.x[2],5:r.x[3]}); res.append(r.fun)
    print(f"seed {sd}: L={r.fun:.4f} r(0,1,4,5)={np.round(r.x,3)} fm; a={o[0]:.3f}; perdas proton {o[2]:.2f} neutron {o[3]:.2f} Rn {o[4]:.2f}",flush=True)
print(f"=> L minimo = {min(res):.4f}; viavel (L<1): {min(res)<1}")
