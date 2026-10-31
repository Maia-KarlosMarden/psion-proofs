import math, itertools, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
ANG=[(-30,30),(90,150),(210,270)]
def r2(t,X):
    Y=np.array([X[k] for k in t]); q=Q9[list(t)]; c=Y.mean(0)
    v=float(np.sum(q*np.sum((Y-c)**2,1))); Q=q.sum(); return v/Q if abs(Q)>1e-12 else v
def build(pairs):
    X={}
    for (a,b),(t1,t2) in zip(pairs,ANG):
        X[a]=np.array([math.cos(math.radians(t1)),math.sin(math.radians(t1))]); X[b]=np.array([math.cos(math.radians(t2)),math.sin(math.radians(t2))])
    return X
def ev(pairs):
    X=build(pairs); rp=r2((0,4,5),X)
    if rp<=0: return None
    return r2((0,1,5),X)/rp, r2((0,1,2),X)/rp, math.sqrt(0.7071/rp)
def ok(r): return r and abs(r[0]/-0.164-1)<0.1 and abs(r[1]-0.86)<0.21
for name,pr in [("hexagono (pares vizinhos)",[(5,0),(1,2),(3,4)]),("PRINCIPAL espelho",[(0,5),(1,4),(2,3)])]:
    r=ev(pr); print(f"{name}: " + ("r_p^2<=0" if r is None else f"R_n={r[0]:+.4f} R_S={r[1]:+.4f} a={r[2]:.3f} fm -> {'ACERTA' if ok(r) else 'nao acerta'}"))
mir=[(0,5),(1,4),(2,3)]; res=[]
for perm in itertools.permutations(mir):
    for fl in itertools.product((0,1),repeat=3):
        pr=[(p[1],p[0]) if f else p for p,f in zip(perm,fl)]; res.append((pr,ev(pr)))
Rns=sorted(set(round(r[0],4) for _,r in res if r)); print("espelho, 48 atribuicoes: R_n distintos =",Rns)
print("passam:",sum(1 for _,r in res if ok(r)),"de 48; r_p^2<=0 em",sum(1 for _,r in res if r is None))
