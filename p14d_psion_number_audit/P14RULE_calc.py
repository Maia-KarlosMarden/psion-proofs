import itertools
TR={"n":(0,1,5),"S0":(0,2,4),"L":(1,2,3),"Spa":(1,3,5),"Spb":(2,3,4)}
LAB={"n":"H","S0":"H","L":"D","Spa":"O","Spb":"D"}
isR=lambda k:k>=3
F=[(f"has{j}",(lambda t,j=j:j in t)) for j in range(6)]
F+=[("has2or3",lambda t:2 in t or 3 in t),("has2and3",lambda t:2 in t and 3 in t),("sum0mod6",lambda t:sum(t)%6==0),
    ("nR>=2",lambda t:sum(map(isR,t))>=2),("nR==1",lambda t:sum(map(isR,t))==1),("hasLR",lambda t:0<sum(map(isR,t))<3),
    ("even>=2",lambda t:sum(k%2==0 for k in t)>=2),("hasmod1",lambda t:any(k%3==1 for k in t))]
CL="HOD"
def fit(train,lab):
    best=None
    for nm,f in F:
        for a in CL:
            for b in CL:
                acc=sum((b if f(TR[h]) else a)==lab[h] for h in train)
                if best is None or acc>best[0]: best=(acc,nm,f,a,b)
    return best
def loo(lab):
    hs=list(TR); sc=0; maj=0; det=[]
    for h in hs:
        tr=[x for x in hs if x!=h]; _,nm,f,a,b=fit(tr,lab); pred=b if f(TR[h]) else a; sc+=pred==lab[h]
        from collections import Counter
        m=Counter(lab[x] for x in tr).most_common(1)[0][0]; maj+=m==lab[h]; det.append((h,nm,a,b,pred,lab[h]))
    return sc,maj,det
if __name__=="__main__":
    sc,maj,det=loo(LAB)
    for d in det: print("held-out %s: stump %s (0->%s,1->%s) predicts %s, true %s"%d)
    print(f"score LOO regra: {sc}/5; maioria: {maj}/5")
    vals=[]
    for perm in set(itertools.permutations([LAB[h] for h in TR])):
        vals.append(loo(dict(zip(TR,perm)))[0])
    import numpy as np; vals=np.array(vals)
    print(f"permutacoes distintas: {len(vals)}; p(score>=obs) = {(vals>=sc).mean():.3f}; media do score nulo = {vals.mean():.2f}")
    print("SUCESSO" if sc>=4 and (vals>=sc).mean()<=0.05 else "SEM SUCESSO")
