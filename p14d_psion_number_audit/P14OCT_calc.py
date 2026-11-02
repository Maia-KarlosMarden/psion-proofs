import math, itertools, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
Q9=np.array([1,-5,-5,10,4,4])/9; NS=8; HOME=(0,1,2,4,5,6)
TH=np.array([2*math.pi*s/NS for s in range(NS)]); XY=np.stack([np.cos(TH),np.sin(TH)],1)
PL=list(itertools.permutations(range(NS),3)); PA=np.array(PL); IDX={p:i for i,p in enumerate(PL)}
TR={"p":(0,4,5),"n":(0,1,5),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def diag_arr(t):
    q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum(); X=XY[PA]
    d=lambda a,b: np.linalg.norm(X[:,a]-X[:,b],axis=1)
    A=-sum(np.cos(TH[PA[:,a]]-TH[HOME[t[a]]]) for a in range(3))
    B=sum(q[a]*q[b]/d(a,b) for a in range(3) for b in range(a+1,3))
    C=-sum(1/d(a,b) for a in range(3) for b in range(a+1,3) if ch[a]!=ch[b]) if len(set(ch))>1 else np.zeros(len(PL))
    c=X.mean(1,keepdims=True); v=(q[None,:]*((X-c)**2).sum(2)).sum(1); R2=v if abs(Qt)<1e-12 else v/Qt
    return A,B,C,R2
def hop_mats():
    r1=[];c1=[];r2=[];c2=[]
    for i,p in enumerate(PL):
        occ=set(p)
        for a in range(3):
            for dlt in (1,-1,2,-2):
                s=(p[a]+dlt)%NS
                if s not in occ:
                    q=list(p); q[a]=s; r1.append(i); c1.append(IDX[tuple(q)])
            s=(p[a]+4)%NS
            if s not in occ:
                q=list(p); q[a]=s; r2.append(i); c2.append(IDX[tuple(q)])
        for a in range(3):
            for b in range(a+1,3):
                q=list(p); q[a],q[b]=q[b],q[a]; r2.append(i); c2.append(IDX[tuple(q)]); r2.append(IDX[tuple(q)]); c2.append(i)
    n=len(PL)
    T=sp.coo_matrix((np.ones(len(r1)),(r1,c1)),shape=(n,n)).tocsr(); T=((T+T.T)>0).astype(float)
    Sx=sp.coo_matrix((np.ones(len(r2)),(r2,c2)),shape=(n,n)).tocsr(); Sx=((Sx+Sx.T)>0).astype(float)
    return T,Sx
T1,S1=hop_mats()
AR={k:diag_arr(t) for k,t in TR.items()}

from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc
H_=list(TR)

from multiprocessing import Pool
from scipy.stats import qmc
H_=list(TR)

def gse(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C
    if t<1e-9 and s<1e-9:
        o=np.sort(dg)[:2]
        if o[1]-o[0]<1e-9:
            i0=np.where(dg<o[0]+1e-9)[0]
            if np.ptp(R2[i0])>1e-6: return [np.nan]*5
            i=i0[0]
        else: i=int(np.argmin(dg))
        return [float(R2[i]),float(dg[i]),float(A[i]),float(B[i]),float(C[i])]
    Hm=sp.diags(dg)-t*T1-s*S1
    w,v=eigsh(Hm.tocsc(),k=2,which='SA'); o=np.argsort(w); w=w[o]; v=v[:,o]
    if w[1]-w[0]<1e-8 and abs((v[:,0]**2*R2).sum()-(v[:,1]**2*R2).sum())>1e-6: return [np.nan]*5
    p2=v[:,0]**2
    return [float((p2*R2).sum()),float(w[0]),float((p2*A).sum()),float((p2*B).sum()),float((p2*C).sum())]
def work(x):
    t,s,k,l=x; return sum([gse(h,t,s,k,l) for h in H_],[])   # 6 hadrons x 5
def Lrad(R):
    L=np.maximum.reduce([np.abs(R[:,0]+0.164)/0.0164,np.abs(R[:,1])/0.15,np.abs(R[:,2]-0.15)/0.10,np.minimum(np.abs(R[:,3]-0.975),np.abs(R[:,4]-0.975))/0.175])
    return np.where(np.isnan(L),50,L)
def unpack(TAB):
    T=TAB.reshape(len(TAB),6,5); rp=T[:,0,0]
    R=np.stack([T[:,i,0]/rp for i in (1,2,3,4,5)],1)   # n,S0,L,Sp1,Sp2
    return R,T
IDX={h:i for i,h in enumerate(H_)}
def mass_ok(M,r1t=0.437,r2t=1.425):
    D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        r1=(M["S0"]-M["L"])/D; a=(M["Sp1"]-M["p"])/D; b=(M["Sp2"]-M["p"])/D
        return (D>0)&(np.abs(r1/r1t-1)<0.3)&((np.abs(a/r2t-1)<0.3)|(np.abs(b/r2t-1)<0.3))


rngD=np.random.default_rng(29)
Dp=rngD.dirichlet([1,1,1],2000); DIRS=np.vstack([Dp,-Dp]); MX=np.abs(DIRS).max(1)
def qmax(Av):
    M={h:DIRS@Av[h] for h in Av}; m=mass_ok(M)
    if not m.any(): return None
    D=(M["L"]-M["n"]); q=D/MX
    return float(q[m].max())
def evalesc(z):
    t,s,k,l=z
    r={h:gse(h,t,s,k,l) for h in H_}
    if any(np.isnan(r[h][0]) for h in H_): return 50.0,None
    rp=r["p"][0]; R=[r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]
    Lr=max(abs(R[0]+0.164)/0.0164,abs(R[1])/0.15,abs(R[2]-0.15)/0.10,min(abs(R[3]-0.975),abs(R[4]-0.975))/0.175)
    q=qmax({h:np.array(r[h][2:5]) for h in H_})
    if q is None: return 40.0,(R,None,Lr)
    return Lr+20*max(0,1-q/0.336),(R,q,Lr)
LO=np.array([0,0,0,0.]); HI=np.array([6,6,20,20.])
def job(z0):
    f=lambda z: evalesc(np.clip(z,LO,HI))[0]
    r=minimize(f,np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':250,'xatol':1e-4,'fatol':1e-5}); z=np.clip(r.x,LO,HI); return float(r.fun),z
def qrow(i):
    return qmax({h:T[i,IDX[h],2:5] for h in H_})

def pstats(x):
    t,s_,k,l=x
    r={h:gse(h,t,s_,k,l) for h in H_}
    if any(np.isnan(r[h][0]) for h in H_): return [np.nan]*7
    rp=r["p"][0]; R=[r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]
    Lr=max(abs(R[0]+0.164)/0.0164,abs(R[1])/0.15,abs(R[2]-0.15)/0.10,min(abs(R[3]-0.975),abs(R[4]-0.975))/0.175)
    q=qmax({h:np.array(r[h][2:5]) for h in H_}); return [Lr,np.nan if q is None else q]+R
def jobL(z0):
    f=lambda z: pstats(np.clip(z,LO,HI))[0] if not np.isnan(pstats(np.clip(z,LO,HI))[0]) else 50.0
    r=minimize(f,np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':300,'xatol':1e-4,'fatol':1e-5}); z=np.clip(r.x,LO,HI); return float(r.fun),z
def job3(z0):
    f=lambda z: evalesc(np.clip(z,LO,HI))[0]
    r=minimize(f,np.clip(z0,LO,HI),method='Nelder-Mead',options={'maxiter':300,'xatol':1e-4,'fatol':1e-5}); z=np.clip(r.x,LO,HI); return float(r.fun),z
if __name__=="__main__":
    rng=np.random.default_rng(31); N=12000
    smp=qmc.LatinHypercube(d=4,seed=31).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: out=np.array(P.map(pstats,list(X),chunksize=40))
    Lr=out[:,0]; q=out[:,1]; R=out[:,2:]
    print(f"tabela: {N} pontos; validos {int(np.sum(~np.isnan(Lr)))}; q definido {int(np.sum(~np.isnan(q)))}")
    print(f"[2] raios isolados: min L_r = {np.nanmin(Lr):.3f}; L_r<1: {int(np.nansum(Lr<1))}; L_r<2: {int(np.nansum(Lr<2))}")
    print(f"    massa isolada: q_max mediana {np.nanmedian(q):.3f}, p90 {np.nanpercentile(q,90):.3f}, max {np.nanmax(q):.3f}; q>=0,336: {int(np.nansum(q>=0.336))}; q>=0,168: {int(np.nansum(q>=0.168))}")
    for thr in (0.336,0.168,0.05,1e-9):
        m=(q>=thr); v=np.where(m,Lr,np.inf); print(f"[3] min L_r com q>={thr:.3g}: {np.nanmin(v):.3f}  (pontos {int(m.sum())})")
    i=int(np.nanargmin(np.where(q>=0.336,Lr,np.inf))) if np.any(q>=0.336) else None
    if i is not None: print(f"    melhor com q>=0,336 (tabela): L_r={Lr[i]:.2f} q={q[i]:.3f} (t,s,k,l)={np.round(X[i],3)} R={np.round(R[i],3)}")
    sc=np.where(np.isnan(q),40.0,Lr+20*np.maximum(0,1-np.nan_to_num(q)/0.336)); sc=np.where(np.isnan(Lr),50.0,sc)
    st3=[X[i] for i in np.argsort(sc)[:20]]; stL=[X[i] for i in np.argsort(np.where(np.isnan(Lr),99,Lr))[:10]]
    with Pool() as P: r3=P.map(job3,st3); rL=P.map(jobL,stL)
    r3.sort(key=lambda a:a[0]); rL.sort(key=lambda a:a[0])
    print(f"[4] raio puro: min L_r apos refino = {rL[0][0]:.3f} ({'VIAVEL' if rL[0][0]<1 else 'inviavel'}); viaveis {sum(1 for r in rL if r[0]<1)}/{len(rL)}")
    f,z=rL[0]; p_=pstats(z); print(f"    melhor raio puro: (t,s,k,l)=({z[0]:.3f},{z[1]:.3f},{z[2]:.2f},{z[3]:.2f}) q={p_[1]:.4f} R={np.round(p_[2:],3)}")
    print(f"[4] conjunto (escala natural): min perda = {r3[0][0]:.3f} ({'VIAVEL' if r3[0][0]<1 else 'inviavel'}); viaveis {sum(1 for r in r3 if r[0]<1)}/{len(r3)}")
    for f,z in r3[:3]:
        p_=pstats(z); print(f"    perda={f:.3f} L_r={p_[0]:.3f} q={p_[1]:.4f} (t,s,k,l)=({z[0]:.3f},{z[1]:.3f},{z[2]:.2f},{z[3]:.2f}) R={np.round(p_[2:],3)}")
