import math, itertools, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
S=np.array([[math.cos(math.pi*i/6),math.sin(math.pi*i/6)] for i in range(12)])
T={"p":(0,4,5),"n":(0,1,5),"Sm":(0,1,2),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def vals(P):  # P: array (N,6) de indices de sitio -> dict de arrays
    out={}
    for k,t in T.items():
        X=S[P[:,list(t)]]; q=Q9[list(t)]; c=X.mean(1,keepdims=True)
        v=(q*((X-c)**2).sum(2)).sum(1); Qt=q.sum(); out[k]=v/Qt if abs(Qt)>1e-12 else v
    rp=out["p"]; R={k:out[k]/np.where(rp>0,rp,np.nan) for k in T}; return R
def crit(R):
    core=(np.abs(R["n"]/-0.164-1)<0.1)&(np.abs(R["Sm"]-0.86)<0.21)
    sp=((R["Sp1"]>=0.8)&(R["Sp1"]<=1.15))|((R["Sp2"]>=0.8)&(R["Sp2"]<=1.15))
    ext=core&(np.abs(R["S0"])<=0.15)&(R["L"]>0)&(R["L"]<=0.25)&sp
    return core,ext
def show(name,P):
    R=vals(np.array([P])); c,e=crit(R)
    print(name,{k:round(float(v[0]),3) for k,v in R.items()},"nucleo",bool(c[0]),"estendido",bool(e[0]))
show("D1 meio-dodecagono",[0,1,2,3,4,5]); show("D2 alternados    ",[0,2,4,6,8,10])
P=np.array(list(itertools.permutations(range(12),6)))
R=vals(P); c,e=crit(R)
print(f"D3: {len(P)} colocacoes; passam nucleo {int(c.sum())}; passam estendido {int(e.sum())}")
if e.sum():
    idx=np.where(e)[0]
    def canon(p):  # ordem por rotacao/reflexao do dodecagono
        best=None
        for r in range(12):
            for s in (1,-1):
                q=tuple((s*x+r)%12 for x in p); best=q if best is None or q<best else best
        return best
    cl={}
    for i in idx: cl.setdefault(canon(tuple(P[i])),i)
    print(f"   classes distintas (a menos de simetria do dodecagono): {len(cl)}")
    for k,i in list(cl.items())[:12]:
        print("   sitios",k,{kk:round(float(R[kk][i]),3) for kk in T})
