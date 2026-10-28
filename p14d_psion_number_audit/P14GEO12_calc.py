import math, itertools, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
S=np.array([[math.cos(math.pi*i/6),math.sin(math.pi*i/6)] for i in range(12)])
NU={"nu_tau (0,1,5)":(0,1,5),"nu_mu (0,2,4)":(0,2,4),"nu_e (1,2,3)":(1,2,3)}
def r2(P,t):
    X=S[P[:,list(t)]]; q=Q9[list(t)]; c=X.mean(1,keepdims=True); return (q*((X-c)**2).sum(2)).sum(1)
def rp_hex():
    X=np.array([[math.cos(math.pi*k/3),math.sin(math.pi*k/3)] for k in (0,4,5)]); q=Q9[[0,4,5]]; c=X.mean(0)
    return float(np.sum(q*np.sum((X-c)**2,1)))
RP=rp_hex()  # proton hexagonal (escala comum)
P1=np.array([[0,1,2,3,4,5]])
print("T1 D1 (theta=30k): "+", ".join(f"{k}: {r2(P1,t)[0]/RP:+.4f}" for k,t in NU.items()))
P=np.array(list(itertools.permutations(range(12),6)))
Z={k:np.abs(r2(P,t)/RP)<1e-9 for k,t in NU.items()}
for k,z in Z.items(): print(f"T2 zeram {k}: {int(z.sum())}/{len(P)}")
allz=Z["nu_tau (0,1,5)"]&Z["nu_mu (0,2,4)"]&Z["nu_e (1,2,3)"]
print(f"T2 zeram os tres: {int(allz.sum())}")
if allz.sum():
    idx=np.where(allz)[0]
    def canon(p):
        return min(tuple((s*x+r)%12 for x in p) for r in range(12) for s in (1,-1))
    cl=sorted({canon(tuple(P[i])) for i in idx}); print(f"T3 classes (a menos de simetria): {len(cl)}"); 
    for c in cl[:20]: print("   ",c)
    aff=[(a,b) for a in range(1,12) for b in range(12) if allz[np.where((P==np.array([(a*k+b)%12 for k in range(6)])).all(1))[0]].any()] if True else []
    print("T4 regras afins s = a k + b que zeram os tres:",aff[:12])
else:
    # mais perto possivel
    m=np.max(np.stack([np.abs(r2(P,t)/RP) for t in NU.values()]),0); i=int(np.argmin(m))
    print(f"T3 nenhuma; melhor max|R| = {m[i]:.4f} em {tuple(P[i])}")
    print("T4 regras afins: nenhuma (nenhuma colocacao zera os tres)")
