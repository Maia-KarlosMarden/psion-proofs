import itertools, math, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
def R2(t,Q):
    th=2*math.pi*np.array(t)/6; q=Q[list(t)]; X=np.stack([np.cos(th),np.sin(th)],1); c=X.mean(0)
    v=(q*((X-c)**2).sum(1)).sum(); Qt=q.sum(); return v if abs(Qt)<1e-12 else v/Qt
ok=tot=0
for perm in itertools.permutations(range(6)):
    Q=Q9[list(perm)]
    try:
        rp=R2((0,4,5),Q); rs=R2((0,1,2),Q)/rp
    except ZeroDivisionError: continue
    if not np.isfinite(rs): continue
    tot+=1; ok+= abs(rs-0.862)<=0.297
print(f"ADENDO (nao registrado antes, rotulado post hoc): em {tot} permutacoes das cargas Q9 sobre os rotulos, fracao com |R_Sigma- - 0.862| <= 1 sigma = {ok/tot:.3f}")
