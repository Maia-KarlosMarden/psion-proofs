import numpy as np
from multiprocessing import Pool
from scipy.stats import qmc
import P14OOS_calc as O
extra,checks,DIRS,POS=O.extra,O.checks,O.DIRS,O.POS
def dm(Av,dirs,v0):
    M={h:dirs@Av[h] for h in Av}; D=M["L"]-M["n"]; u=176.1/D
    return ((M["Sm"]-M["S0"])*u)
def rowb(z):
    o=extra(z)
    if o is None: return (0,0)
    R,Av,Rsm=o; v0,_,_,_,_=checks(Av,Rsm,POS); d=dm(Av,POS,v0)[v0]
    return int(v0.sum()),int(((d>=0)&(d<=50)).sum())
if __name__=="__main__":
    OPT={"O1":[0.258,0.224,0.047,0.229],"O2":[0.080,0.439,0.0,0.083],"O3":[0.035,0.456,0.095,0.0]}
    for nm,z in OPT.items():
        R,Av,Rsm=extra(np.array(z)); v0,_,_,_,_=checks(Av,Rsm,DIRS); d=dm(Av,DIRS,v0)[v0]
        ok=(d>=0)&(d<=50)
        print(f"{nm}: V0 dirs {int(v0.sum())}; dM range [{d.min():.0f},{d.max():.0f}] MeV (obs 4.8); P-C pass {int(ok.sum())} ({ok.mean():.2f})")
    rng=np.random.default_rng(53); N=1500; smp=qmc.LatinHypercube(d=4,seed=53).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: res=np.array(P.map(rowb,list(X),chunksize=30))
    print(f"baseline: V0 pairs {int(res[:,0].sum())}; P-C pass fraction p0 = {res[:,1].sum()/max(res[:,0].sum(),1):.3f}")
