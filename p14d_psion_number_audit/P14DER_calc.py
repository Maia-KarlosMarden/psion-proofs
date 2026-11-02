import numpy as np, math
import P14FF3_calc as F
hc=0.19733
def Lcommon(r):
    o=F.evalr({l:r for l in range(6)}); return None if o is None else (max(o[2:]),o)
print("C1 (r=1.6e-28 fm ~ 0):",F.Lof({l:0.0 for l in range(6)}))
r2=hc/0.610; L=Lcommon(r2); print(f"C2 r={r2:.4f} fm:", None if L is None else round(L[0],3))
for nm,r in (("C3a",hc/0.341),("C3b",2.04*hc/0.341)):
    L=Lcommon(r); print(f"{nm} r={r:.3f} fm:", "a^2<0 (invalido)" if L is None else round(L[0],3))
best=(1e9,None)
for c in np.linspace(0.1,2.0,191):
    m={0:c/2.5,1:c/1.5,4:c/1.5,5:c/2.5,2:c/0.5,3:c/0.5}
    o=F.evalr(m)
    if o is None: continue
    L=max(o[2:]); best=min(best,(L,c))
c=best[1]; m={0:c/2.5,1:c/1.5,4:c/1.5,5:c/2.5}; o=F.evalr(m)
print(f"C4 melhor c = {c:.2f} fm: L = {best[0]:.3f}; tamanhos r(0,1,4,5) = {c/2.5:.3f},{c/1.5:.3f},{c/1.5:.3f},{c/2.5:.3f} fm; a={o[0]:.3f}; perdas p {o[2]:.2f} n {o[3]:.2f} Rn {o[4]:.2f}")
