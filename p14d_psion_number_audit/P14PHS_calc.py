import math, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
TH=[math.pi*k/3 for k in range(6)]
def wrap(x): return (x+math.pi)%(2*math.pi)-math.pi
def frechet(t):
    best=None
    for m in np.linspace(-math.pi,math.pi,36001):
        f=sum(wrap(TH[k]-m)**2 for k in t)
        if best is None or f<best[0]-1e-12: best=(f,m)
    ms=[m for m in np.linspace(-math.pi,math.pi,36001) if abs(sum(wrap(TH[k]-m)**2 for k in t)-best[0])<1e-9]
    # minimo unico se os minimizadores formam um unico aglomerado
    spread=max(abs(wrap(x-ms[0])) for x in ms); return best[1], len(ms)>0 and spread<0.05   # distancia circular: corrige duplicacao em m=+-pi
def mom(t,met):
    q=Q9[list(t)]; Q=q.sum()
    if met=="I1":
        X=np.array([[math.cos(TH[k]),math.sin(TH[k])] for k in t]); c=X.mean(0); v=float(np.sum(q*np.sum((X-c)**2,1)))
    elif met=="I2":
        m,ok=frechet(t)
        if not ok: return None
        v=float(sum(q[i]*wrap(TH[k]-m)**2 for i,k in enumerate(t)))
    else:
        z=sum(np.exp(1j*TH[k]) for k in t)
        if abs(z)<1e-9: return None
        m=np.angle(z); v=float(sum(q[i]*abs(np.exp(1j*TH[k])-np.exp(1j*m))**2 for i,k in enumerate(t)))
    return v/Q if abs(Q)>1e-12 else v
T={"p":(0,4,5),"n":(0,1,5),"Sigma-":(0,1,2),"Sigma0":(0,2,4),"Lambda":(1,2,3),"Sigma+ (1,3,5)":(1,3,5),"Sigma+ (2,3,4)":(2,3,4),"(3,4,5)":(3,4,5)}
ok={"n":lambda R:abs(R/-0.164-1)<0.1,"Sigma-":lambda R:abs(R-0.86)<0.21,"Sigma0":lambda R:abs(R)<=0.15,"Lambda":lambda R:0<R<=0.25,"Sigma+ (1,3,5)":lambda R:0.8<=R<=1.15,"Sigma+ (2,3,4)":lambda R:0.8<=R<=1.15}
for met,nm in (("I1","I1 corda/centroide (controle = hexagono)"),("I2","I2 arco/Frechet"),("I3","I3 corda/media circular")):
    rp=mom(T["p"],met); print(nm)
    npass=0
    for k,t in T.items():
        if k=="p": continue
        v=mom(t,met)
        if v is None: print(f"   {k:16s}: indefinida (sem previsao)"); continue
        R=v/rp; res=ok.get(k); s=("" if res is None else (" passa" if res(R) else " falha")); npass+= (res is not None and res(R))
        print(f"   {k:16s}: R = {R:+.3f}{s}")
    print(f"   passam: {npass}")
