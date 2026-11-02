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

from multiprocessing import Pool
from scipy.stats import qmc
H_=list(TR)
def gse(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C
    if t<1e-9 and s<1e-9:
        o=np.sort(dg)[:2]
        if o[1]-o[0]<1e-9:
            i0=np.where(dg<o[0]+1e-9)[0]
            return (float(R2[i0[0]]),float(o[0])) if np.ptp(R2[i0])<1e-6 else (np.nan,np.nan)
        return float(R2[np.argmin(dg)]),float(o[0])
    Hm=sp.diags(dg)-t*T1-s*S1
    w,v=eigsh(Hm.tocsc(),k=2,which='SA'); o=np.argsort(w); w=w[o]; v=v[:,o]
    if w[1]-w[0]<1e-8 and abs((v[:,0]**2*R2).sum()-(v[:,1]**2*R2).sum())>1e-6: return np.nan,np.nan
    return float((v[:,0]**2*R2).sum()),float(w[0])
def work(x):
    t,s,k,l=x; r={h:gse(h,t,s,k,l) for h in H_}; rp=r["p"][0]
    return [r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]+[r[h][1] for h in H_]

def analyze(X,TAB):
    R=TAB[:,:5]; E={h:TAB[:,5+i] for i,h in enumerate(H_)}
    D=E["L"]-E["n"]
    with np.errstate(all='ignore'):
        r1=(E["S0"]-E["L"])/D; sp1=(E["Sp1"]-E["p"])/D; sp2=(E["Sp2"]-E["p"])/D
        r1ok=np.abs(r1/0.437-1)<0.3; r2ok=(np.abs(sp1/1.425-1)<0.3)|(np.abs(sp2/1.425-1)<0.3)
        c=1.29/(176.1/D)-(E["n"]-E["p"]); cD=np.abs(c)/D
    ok=(D>0)&r1ok&r2ok
    L=np.maximum.reduce([np.abs(R[:,0]+0.164)/0.0164,np.abs(R[:,1])/0.15,np.abs(R[:,2]-0.15)/0.10,np.minimum(np.abs(R[:,3]-0.975),np.abs(R[:,4]-0.975))/0.175])
    L=np.where(np.isnan(L),50,L)
    print(f"pontos {len(X)}: D>0: {int((D>0).sum())}; R1: {int((D>0)&r1ok).sum() if False else int(((D>0)&r1ok).sum())}; R2: {int(((D>0)&r2ok).sum())}; D>0&R1&R2: {int(ok.sum())} ({100*ok.mean():.2f}%)")
    if ok.sum():
        i=np.where(ok)[0][np.argmin(L[ok])]
        print(f"   entre eles: melhor L de raios = {L[ok].min():.2f}; L<2: {int((L[ok]<2).sum())}; L<1: {int((L[ok]<1).sum())}")
        print(f"   melhor em (t,s,k,l)=({X[i][0]:.3f},{X[i][1]:.3f},{X[i][2]:.2f},{X[i][3]:.2f}) R=" + ",".join(f"{v:+.3f}" for v in R[i]) + f"; |c|/D={cD[i]:.2f}; c u = {c[i]*176.1/D[i]:+.1f} MeV")
        print(f"   |c|/D entre os que passam: mediana {np.median(cD[ok]):.2f}; com |c|/D<=3: {int((ok&(cD<=3)).sum())}")
        j=ok&(L<2)
        if j.sum():
            for i in np.where(j)[0][:5]: print(f"   L<2: L={L[i]:.2f} (t,s,k,l)=({X[i][0]:.3f},{X[i][1]:.3f},{X[i][2]:.2f},{X[i][3]:.2f}) |c|/D={cD[i]:.2f}")
    print(f"   (controle) L<1 sem massa: {int((L<1).sum())}; min L = {L.min():.2f}")
if __name__=="__main__":
    rng=np.random.default_rng(25); N=12000
    smp=qmc.LatinHypercube(d=4,seed=25).random(N)
    zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: TAB=np.array(P.map(work,list(X),chunksize=40))
    np.savez('P14MASS_table.npz',X=X,TAB=TAB)
    analyze(X,TAB)
