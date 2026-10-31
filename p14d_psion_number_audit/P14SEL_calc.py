import math, itertools, numpy as np
rng=np.random.default_rng(18)
Q9=np.array([1,-5,-5,10,4,4])/9
TH=np.array([math.pi/3*s for s in range(6)]); XY=np.stack([np.cos(TH),np.sin(TH)],1)
PA=np.array(list(itertools.permutations(range(6),3)))
TR={"p":(0,4,5),"n":(0,1,5),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def arr(t):
    q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum(); X=XY[PA]
    d=lambda a,b: np.linalg.norm(X[:,a]-X[:,b],axis=1)
    A=-sum(np.cos(TH[PA[:,a]]-TH[t[a]]) for a in range(3))
    B=sum(q[a]*q[b]/d(a,b) for a in range(3) for b in range(a+1,3))
    C=-sum(1/d(a,b) for a in range(3) for b in range(a+1,3) if ch[a]!=ch[b]) if len(set(ch))>1 else np.zeros(len(PA))
    c=X.mean(1,keepdims=True); v=(q[None,:]*((X-c)**2).sum(2)).sum(1); R2=v if abs(Qt)<1e-12 else v/Qt
    return A,B,C,R2
AR={k:arr(t) for k,t in TR.items()}
ks=np.arange(0,30.0001,1.0); GK,GL=np.meshgrid(ks,ks,indexing='ij'); GK=GK.ravel(); GL=GL.ravel(); G=len(GK)
Ts=[0.1,0.3,1,3,10,np.inf]
def avg(k,T):
    A,B,C,R2=AR[k]; E=A[None,:]+GK[:,None]*B[None,:]+GL[:,None]*C[None,:]
    if np.isinf(T): w=np.ones_like(E)
    else: w=np.exp(-(E-E.min(1,keepdims=True))/T)
    return (w*R2[None,:]).sum(1)/w.sum(1)
W=(0.0164,0.15,0.10,0.175)
RR={}
for T in Ts:
    rp=avg("p",T); RR[T]={k:avg(k,T)/rp for k in ("n","S0","L","Sp1","Sp2")}
def feas(T,cn,cs0,cl,csp):
    R=RR[T]; return (np.abs(R["n"]-cn)<W[0])&(np.abs(R["S0"]-cs0)<W[1])&(np.abs(R["L"]-cl)<W[2])&((np.abs(R["Sp1"]-csp)<W[3])|(np.abs(R["Sp2"]-csp)<W[3]))
tot=0
for T in Ts:
    R=RR[T]; f=feas(T,-0.164,0.0,0.15,0.975); tot+=int(f.sum())
    ind={nm:int(m.sum()) for nm,m in (("n",np.abs(R["n"]+0.164)<W[0]),("S0",np.abs(R["S0"])<W[1]),("L",np.abs(R["L"]-0.15)<W[2]),("Sp",(np.abs(R["Sp1"]-0.975)<W[3])|(np.abs(R["Sp2"]-0.975)<W[3])))}
    print(f"T={T}: completos {int(f.sum())}/{G}; individuais {ind}")
print(f"raios completos total: {tot} de {G*len(Ts)}")
hits=0;N=2000
for _ in range(N):
    c=rng.uniform(-0.6,2.2,4)
    if any(feas(T,*c).any() for T in Ts): hits+=1
print(f"teste de acaso: P_acaso = {hits/N:.3f} ({hits}/{N})")
