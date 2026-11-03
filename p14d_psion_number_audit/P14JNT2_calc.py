import math, numpy as np
from scipy.optimize import differential_evolution
import P14SPINF_calc as SF
import P14FF3_calc as F
hc=0.19733
SPIN={}
for rule in ("R1","R4"):
    SPIN[rule]={nm:SF.state(t,rule) for nm,t in (("p",(0,4,5)),("n",(0,1,5)))}
def terms(x,rule):
    sz={0:x[0],1:x[1],4:x[2],5:x[3]}
    o=F.evalr(sz)
    if o is None: return None
    a=o[0]; te=[o[2],o[3],o[4]]
    out=[]
    for nm,t in (("p",(0,4,5)),("n",(0,1,5))):
        lab,qq,d=F.geo(t); s=SPIN[rule][nm]
        G=lambda q: sum(qi*si*math.exp(-q*q*sz[l]**2/6)*F.j0(q*a*di) for l,qi,di,si in zip(lab,qq,d,s))
        g0=G(0.0)
        if abs(g0)<1e-9: return None
        out.append(max(abs(G(math.sqrt(Q2)/hc)/g0/(1+Q2/0.71)**-2-1) for Q2 in F.Q2P)/0.15)
    mu=lambda nm,t: sum(F.Q9[l]*si for l,si in zip(t,SPIN[rule][nm]))
    z6=abs(mu("n",(0,1,5))/mu("p",(0,4,5))+0.685)/0.20
    return te+out+[z6]
def loss(x,rule):
    t=terms(x,rule); return 50.0 if t is None else max(t)
names=["E_p forma","E_n forma","R_n","M_p forma","M_n forma","mu_n/mu_p"]
if __name__=="__main__":
    for rule in ("R1","R4"):
        best=(1e9,None)
        for sd in range(91,101):
            r=differential_evolution(loss,[(0,1.5)]*4,args=(rule,),popsize=40,maxiter=300,seed=sd,tol=1e-10,updating='immediate',polish=True)
            if r.fun<best[0]: best=(r.fun,r.x)
        t=terms(best[1],rule)
        print(f"{rule}: melhor L = {best[0]:.3f}; r(0,1,4,5) = {np.round(best[1],3)} fm; termos: "+", ".join(f"{n} {v:.2f}" for n,v in zip(names,t))+f"; viavel = {best[0]<1}",flush=True)
