import math, itertools, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from scipy.stats import qmc
Q9=np.array([1,-5,-5,10,4,4])/9
PL=list(itertools.permutations(range(6),3)); PA=np.array(PL); NP=len(PL); IDX={p:i for i,p in enumerate(PL)}
TR={"p":(0,4,5),"n":(0,1,5),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}; H_=list(TR)
PHI=np.array([2*math.pi/3*(s%3) for s in range(6)]); LAY=np.array([s//3 for s in range(6)])
# hopping (independente de h)
rt=[];ct=[];rs=[];cs=[]
for i,p in enumerate(PL):
    occ=set(p)
    for a in range(3):
        for dd in (1,2):
            s_=(p[a]//3)*3+((p[a]%3)+dd)%3
            if s_ not in occ: q=list(p); q[a]=s_; rt.append(i); ct.append(IDX[tuple(q)])
        s_=(p[a]+3)%6
        if s_ not in occ: q=list(p); q[a]=s_; rs.append(i); cs.append(IDX[tuple(q)])
    for a in range(3):
        for b in range(a+1,3):
            if p[a]%3==p[b]%3 and p[a]//3!=p[b]//3:
                q=list(p); q[a],q[b]=q[b],q[a]; rs.append(i); cs.append(IDX[tuple(q)])
def mat(r,c):
    M=np.zeros((NP,NP)); M[r,c]=1; return ((M+M.T)>0).astype(float)
T1=mat(rt,ct); S1=mat(rs,cs)
def arrays(t,h):
    q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum()
    xyz=np.stack([np.cos(PHI),np.sin(PHI),h*LAY],1); X=xyz[PA]   # (NP,3,3)
    d=lambda a,b: np.linalg.norm(X[:,a]-X[:,b],axis=1)
    A=-sum((np.cos(PHI[PA[:,a]]-PHI[t[a]])+np.where(LAY[PA[:,a]]==LAY[t[a]],1.0,-1.0))/2 for a in range(3))
    B=sum(q[a]*q[b]/d(a,b) for a in range(3) for b in range(a+1,3))
    C=-sum(1/d(a,b) for a in range(3) for b in range(a+1,3) if ch[a]!=ch[b]) if len(set(ch))>1 else np.zeros(NP)
    c=X.mean(1,keepdims=True); v=(q[None,:]*((X-c)**2).sum(2)).sum(1); R2=v if abs(Qt)<1e-12 else v/Qt
    return A,B,C,np.where(np.abs(R2)<1e-12,0.0,R2)
def gse(h_,t,s,k,l,AR):
    A,B,C,R2=AR[h_]; dg=A+k*B+l*C
    H=np.diag(dg)-t*T1-s*S1
    w,v=np.linalg.eigh(H); p2=v[:,0]**2
    if w[1]-w[0]<1e-8 and abs((v[:,0]**2*R2).sum()-(v[:,1]**2*R2).sum())>1e-6: return [np.nan]*5
    return [float((p2*R2).sum()),float(w[0]),float((p2*A).sum()),float((p2*B).sum()),float((p2*C).sum())]
rngD=np.random.default_rng(33); Dp=rngD.dirichlet([1,1,1],2000); DIRS=np.vstack([Dp,-Dp]); MX=np.abs(DIRS).max(1)
def mass_ok(M,r1t=0.437,r2t=1.425):
    D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        r1=(M["S0"]-M["L"])/D; a=(M["Sp1"]-M["p"])/D; b=(M["Sp2"]-M["p"])/D
        return (D>0)&(np.abs(r1/r1t-1)<0.3)&((np.abs(a/r2t-1)<0.3)|(np.abs(b/r2t-1)<0.3))
def qmax(Av):
    M={h:DIRS@Av[h] for h in Av}; m=mass_ok(M)
    return None if not m.any() else float(((M["L"]-M["n"])/MX)[m].max())
def pst(x):
    t,s,k,l,h=x; AR={n:arrays(tr,h) for n,tr in TR.items()}
    r={n:gse(n,t,s,k,l,AR) for n in H_}
    if any(np.isnan(r[n][0]) for n in H_): return [np.nan]*7
    rp=r["p"][0]; R=[r[n][0]/rp for n in ("n","S0","L","Sp1","Sp2")]
    Lr=max(abs(R[0]+0.164)/0.0164,abs(R[1])/0.15,abs(R[2]-0.15)/0.10,min(abs(R[3]-0.975),abs(R[4]-0.975))/0.175)
    q=qmax({n:np.array(r[n][2:5]) for n in H_}); return [Lr,np.nan if q is None else q]+R
LO=np.array([0,0,0,0,0.4]); HI=np.array([6,6,20,20,3.])
def lossJ(z):
    p=pst(np.clip(z,LO,HI))
    if np.isnan(p[0]): return 50.0
    if np.isnan(p[1]): return 40.0
    return p[0]+20*max(0,1-p[1]/0.336)
def lossR(z):
    p=pst(np.clip(z,LO,HI)); return 50.0 if np.isnan(p[0]) else p[0]
def nm(arg):
    kind,z0=arg; f=lossJ if kind=="J" else lossR
    r=minimize(f,np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':400,'xatol':1e-4,'fatol':1e-5}); return kind,float(r.fun),np.clip(r.x,LO,HI)
if __name__=="__main__":
    rng=np.random.default_rng(34); N=30000
    smp=qmc.LatinHypercube(d=5,seed=34).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u,a,b: np.exp(np.log(a)+u*(np.log(b)-np.log(a)))
    X=np.stack([np.where(zt,0,lg(smp[:,0],0.01,6)),np.where(zs,0,lg(smp[:,1],0.01,6)),smp[:,2]*20,smp[:,3]*20,lg(smp[:,4],0.4,3)],1)
    with Pool() as P: out=np.array(P.map(pst,list(X),chunksize=100))
    Lr=out[:,0]; q=out[:,1]; R=out[:,2:]
    print(f"tabela: {N} pontos; validos {int(np.sum(~np.isnan(Lr)))}; q definido {int(np.sum(~np.isnan(q)))}")
    print(f"[2] raios isolados: min L_r = {np.nanmin(Lr):.3f}; L_r<1: {int(np.nansum(Lr<1))}; L_r<2: {int(np.nansum(Lr<2))}")
    print(f"    massa isolada: q_max mediana {np.nanmedian(q):.3f}, p90 {np.nanpercentile(q,90):.3f}, max {np.nanmax(q):.3f}; q>=0,336: {int(np.nansum(q>=0.336))}; q>=0,168: {int(np.nansum(q>=0.168))}")
    for thr in (0.336,0.168,0.05,1e-9):
        m=(q>=thr); v=np.where(m,Lr,np.inf); print(f"[3] min L_r com q>={thr:.3g}: {np.nanmin(v):.3f} (pontos {int(m.sum())})")
    sc=np.where(np.isnan(Lr),50.0,np.where(np.isnan(q),40.0,Lr+20*np.maximum(0,1-np.nan_to_num(q)/0.336)))
    stJ=[("J",X[i]) for i in np.argsort(sc)[:20]]; stR=[("R",X[i]) for i in np.argsort(np.where(np.isnan(Lr),99,Lr))[:10]]
    with Pool() as P: res=P.map(nm,stJ+stR)
    rJ=sorted([r for r in res if r[0]=="J"],key=lambda a:a[1]); rR=sorted([r for r in res if r[0]=="R"],key=lambda a:a[1])
    print(f"[4] raio puro: min L_r = {rR[0][1]:.3f} ({'VIAVEL' if rR[0][1]<1 else 'inviavel'}); viaveis {sum(1 for r in rR if r[1]<1)}/{len(rR)}")
    z=rR[0][2]; p_=pst(z); print(f"    melhor: (t,s,k,l,h)={np.round(z,3)} q={p_[1]:.4f} R={np.round(p_[2:],3)}")
    print(f"[4] conjunto (escala natural): min perda = {rJ[0][1]:.3f} ({'VIAVEL' if rJ[0][1]<1 else 'inviavel'}); viaveis {sum(1 for r in rJ if r[1]<1)}/{len(rJ)}")
    for _,f,z in rJ[:3]:
        p_=pst(z); print(f"    perda={f:.3f} L_r={p_[0]:.3f} q={p_[1]:.4f} (t,s,k,l,h)={np.round(z,3)} R={np.round(p_[2:],3)}")
