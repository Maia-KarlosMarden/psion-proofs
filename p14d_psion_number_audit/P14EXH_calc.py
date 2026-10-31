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
if __name__=="__main__":
    rng=np.random.default_rng(21)
    smp=qmc.LatinHypercube(d=4,seed=21).random(6000)
    def logu(u,zero):
        return np.where(zero,0.0,np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01))))
    zt=rng.random(6000)<0.05; zs=rng.random(6000)<0.05
    X=np.stack([logu(smp[:,0],zt),logu(smp[:,1],zs),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: Ls=np.array(P.map(loss,list(X),chunksize=30))
    o=np.argsort(Ls); print(f"(1) LHS 6000 pontos: min L = {Ls.min():.3f}; viaveis (L<1): {int((Ls<1).sum())}; L<1.5: {int((Ls<1.5).sum())}")
    for i in o[:5]: show("LHS",X[i],Ls[i])
    best=[]
    for sd in (1,2,3):
        r=differential_evolution(loss,[(0,6),(0,6),(0,20),(0,20)],seed=sd,popsize=12,maxiter=60,tol=0,workers=2,updating='deferred',polish=False)
        print(f"(2) DE semente {sd}: L={r.fun:.3f}"); show("DE",r.x,r.fun); best.append((r.fun,r.x))
    starts=[X[i] for i in o[:8]]+[b[1] for b in best]
    # pontos n+Lambda da grade SWP
    tv=[0,0.03,0.06,0.1,0.15,0.2,0.3,0.5,0.8,1.2,2,3.5]; kv=[0,0.5,1,2,3,5,8,12]
    cand=[]
    for t in tv:
        for s_ in tv:
            for k in kv:
                for l in kv: cand.append((t,s_,k,l))
    with Pool() as P: cl=P.map(loss,cand,chunksize=40)
    sw=[cand[i] for i in np.argsort(cl)[:4]]
    print("(3) 4 melhores da grade SWP (por L):"); 
    for c in sw: show("SWP",np.array(c),loss(c))
    starts+= [np.array(c) for c in sw]
    def nm(x0):
        r=minimize(lambda z: loss(np.clip(z,[0,0,0,0],[6,6,20,20])),x0,method='Nelder-Mead',options={'maxiter':200,'xatol':1e-4,'fatol':1e-5}); z=np.clip(r.x,[0,0,0,0],[6,6,20,20]); return loss(z),z
    with Pool() as P: res=P.map(nm,starts)
    res.sort(key=lambda a:a[0]); print(f"(3) Nelder-Mead de {len(starts)} partidas: min L = {res[0][0]:.3f}")
    for L,z in res[:5]: show("NM",z,L)
    Lmin=min(Ls.min(),min(b[0] for b in best),res[0][0]); print(f"MIN L GLOBAL = {Lmin:.3f}  ->  "+("VIAVEL ENCONTRADO" if Lmin<1 else "esgotado: nenhum ponto viavel em [0,6]^2x[0,20]^2"))
