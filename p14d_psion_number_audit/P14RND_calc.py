import math, itertools, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
Q9=np.array([1,-5,-5,10,4,4])/9; NS=12
TH=np.array([2*math.pi*s/NS for s in range(NS)]); XY=np.stack([np.cos(TH),np.sin(TH)],1)
PL=list(itertools.permutations(range(NS),3)); PA=np.array(PL); IDX={p:i for i,p in enumerate(PL)}
TR={"p":(0,4,5),"n":(0,1,5),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def diag_arr(t):
    q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum(); X=XY[PA]
    d=lambda a,b: np.linalg.norm(X[:,a]-X[:,b],axis=1)
    A=-sum(np.cos(TH[PA[:,a]]-TH[2*t[a]]) for a in range(3))
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
            s=(p[a]+6)%NS
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
def gsr(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C
    if t<1e-9 and s<1e-9:
        o=np.sort(dg)[:2]
        if o[1]-o[0]<1e-9:
            i0=np.where(dg<o[0]+1e-9)[0]
            return float(R2[i0[0]]) if np.ptp(R2[i0])<1e-6 else np.nan
        return float(R2[np.argmin(dg)])
    Hm=sp.diags(dg)-t*T1-s*S1
    w,v=eigsh(Hm.tocsc(),k=2,which='SA'); o=np.argsort(w); w=w[o]; v=v[:,o]
    if w[1]-w[0]<1e-8 and abs((v[:,0]**2*R2).sum()-(v[:,1]**2*R2).sum())>1e-6: return np.nan
    return float((v[:,0]**2*R2).sum())
def Rs(x):
    t,s,k,l=x; r={h:gsr(h,t,s,k,l) for h in H_}; rp=r["p"]
    return {h:r[h]/rp for h in H_ if h!="p"}
def loss(x):
    R=Rs(x)
    if any(np.isnan(v) for v in R.values()): return 50.0
    return max(abs(R["n"]+0.164)/0.0164,abs(R["S0"])/0.15,abs(R["L"]-0.15)/0.10,min(abs(R["Sp1"]-0.975),abs(R["Sp2"]-0.975))/0.175)
def show(tag,x,L):
    R=Rs(x); print(f"   {tag}: L={L:.3f} (t,s,kappa,lambda)=({x[0]:.3f},{x[1]:.3f},{x[2]:.2f},{x[3]:.2f}) "+" ".join(f"{h}={R[h]:+.3f}" for h in R))

def Rv(x):
    R=Rs(x); return [R["n"],R["S0"],R["L"],R["Sp1"],R["Sp2"]]
def Lt(R,c):
    return max(abs(R[0]-c[0])/0.0164,abs(R[1]-c[1])/0.15,abs(R[2]-c[2])/0.10,min(abs(R[3]-c[3]),abs(R[4]-c[3]))/0.175)
LO=np.array([0,0,0,0.]); HI=np.array([6,6,20,20.])
def refine(arg):
    x0,c=arg
    def f(z):
        z=np.clip(z,LO,HI); R=Rv(z)
        return 50.0 if any(np.isnan(R)) else Lt(R,c)
    r=minimize(f,x0,method='Nelder-Mead',options={'maxiter':60,'xatol':1e-3,'fatol':1e-4}); return r.fun
def run_target(TAB,X,c):
    L0=np.array([50.0 if np.any(np.isnan(r)) else Lt(r,c) for r in TAB]); idx=np.argsort(L0)[:2]
    return min(min(L0[idx]),*[refine((X[i],c)) for i in idx])
if __name__=="__main__":
    rng=np.random.default_rng(23); N=12000
    smp=qmc.LatinHypercube(d=4,seed=23).random(N)
    zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: TAB=np.array(P.map(Rv,list(X),chunksize=40))
    print(f"tabela: {N} pontos; NaN: {int(np.isnan(TAB).any(1).sum())}")
    real=(-0.164,0.0,0.15,0.975)
    L0=np.array([50.0 if np.any(np.isnan(r)) else Lt(r,real) for r in TAB]); print(f"alvo real: melhor L da tabela = {L0.min():.3f}")
    Lr=run_target(TAB,X,real); print(f"alvo real: min L apos refino = {Lr:.3f} -> {'SUCESSO' if Lr<1 else 'falha'}")
    M=80; cs=[tuple(rng.uniform(-0.6,2.2,4)) for _ in range(M)]
    def job(c): return run_target(TAB,X,c)
    with Pool() as P: Ls=P.map(job,cs)
    Ls=np.array(Ls); k=int((Ls<1).sum()); p=k/M
    se=1.96*np.sqrt(p*(1-p)/M)
    print(f"alvos aleatorios: {k}/{M} com sucesso -> P_acaso = {p:.3f} (IC95% aprox. [{max(0,p-se):.3f}, {min(1,p+se):.3f}]); mediana de min L = {np.median(Ls):.2f}")
