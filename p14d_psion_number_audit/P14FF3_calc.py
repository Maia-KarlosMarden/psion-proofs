import math, numpy as np
from scipy.optimize import differential_evolution
Q9=np.array([1,-5,-5,10,4,4])/9; hc=0.19733; mN=0.9389; mun=-1.913
def X(t): th=2*math.pi*np.array(t)/6; return np.stack([np.cos(th),np.sin(th)],1)
def geo(t): x=X(t); c=x.mean(0); return list(t),Q9[list(t)],np.linalg.norm(x-c,axis=1)
GP,GN=geo((0,4,5)),geo((0,1,5))
j0=lambda z: 1.0 if abs(z)<1e-8 else math.sin(z)/z
def Fh(g,rmap,a,q): lab,qq,d=g; return sum(qi*math.exp(-q*q*rmap[l]**2/6)*j0(q*a*di) for l,qi,di in zip(lab,qq,d))
def slope_num(g,rmap,a): lab,qq,d=g; return sum(qi*(a*a*di*di+rmap[l]**2) for l,qi,di in zip(lab,qq,d))
Q2P=(0.1,0.2,0.3,0.5,0.8,1.0)
def evalr(rmap):
    lab,qq,d=GP; ring_p=sum(qi*di*di for qi,di in zip(qq,d))/qq.sum(); extra_p=sum(qi*rmap[l]**2 for l,qi in zip(lab,qq))/qq.sum()
    a2=(0.8414**2-extra_p)/ring_p
    if a2<=0: return None
    a=math.sqrt(a2); Rn=slope_num(GN,rmap,a)/0.8414**2
    rows=[]
    for Q2 in Q2P:
        q=math.sqrt(Q2)/hc; GD=(1+Q2/0.71)**-2; tau=Q2/(4*mN**2); GEn=-mun*tau/(1+5.6*tau)*GD
        rows.append((Fh(GP,rmap,a,q)/GD, Fh(GN,rmap,a,q),GEn))
    mx=max(abs(r[2]) for r in rows)
    return a,Rn,max(abs(r[0]-1) for r in rows)/0.15,max(abs(r[1]-r[2]) for r in rows)/(0.15*mx),abs(Rn+0.164)/0.0164
def Lof(rmap):
    o=evalr(rmap)
    return 50.0 if o is None else max(o[2],o[3],o[4])
if __name__=="__main__":
    print("A: tamanho comum")
    best=(1e9,None)
    for r in np.linspace(0,0.9,91):
        o=evalr({l:r for l in range(6)})
        if o is None: continue
        L=max(o[2:]); best=min(best,(L,r))
        if abs(r*100%10)<1e-6 or r in (0,): print(f"  r_psi={r:.2f} fm: a={o[0]:.3f} R_n={o[1]:+.3f}; perda proton {o[2]:.2f}, neutron {o[3]:.2f}, R_n {o[4]:.2f}; L={L:.2f}")
    print(f"  melhor L = {best[0]:.3f} em r_psi = {best[1]:.2f} fm; veredito A: {'SUPPORT' if best[0]<1 else 'KILL'}")
    print("B: tamanhos por psion (r0,r1,r4,r5)")
    f=lambda x: Lof({0:x[0],1:x[1],4:x[2],5:x[3]})
    for sd in range(72,77):
        r=differential_evolution(f,[(0,0.9)]*4,popsize=15,maxiter=40,seed=sd,tol=1e-8,updating='immediate',polish=False)
        o=evalr({0:r.x[0],1:r.x[1],4:r.x[2],5:r.x[3]}); big,small=r.x.max(),r.x[r.x>0.02].min() if (r.x>0.02).any() else float('nan')
        print(f"  seed {sd}: L={r.fun:.3f} r(0,1,4,5)={np.round(r.x,3)} fm; a={o[0]:.3f}; max/min(>0.02)={big/small:.1f}; viavel={r.fun<1}")
