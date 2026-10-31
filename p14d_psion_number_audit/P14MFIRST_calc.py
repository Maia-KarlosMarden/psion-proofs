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
def analyze(X,TAB,tag):
    R=TAB[:,:5]; E={h:TAB[:,5+i] for i,h in enumerate(H_)}
    dn=E["n"]-E["p"]; dL=E["L"]-E["p"]
    with np.errstate(all='ignore'):
        u=177.4/dL; sp1=u*(E["Sp1"]-E["p"]); sp2=u*(E["Sp2"]-E["p"]); s0=u*(E["S0"]-E["L"]); nm=u*dn
    m1=dn>0; m2=dL>0
    m3=(np.abs(sp1/251-1)<0.3)|(np.abs(sp2/251-1)<0.3); m4=np.abs(s0/77-1)<0.3
    m5=(nm>0.65)&(nm<2.6)
    L=np.maximum.reduce([np.abs(R[:,0]+0.164)/0.0164,np.abs(R[:,1])/0.15,np.abs(R[:,2]-0.15)/0.10,np.minimum(np.abs(R[:,3]-0.975),np.abs(R[:,4]-0.975))/0.175])
    L=np.where(np.isnan(L),50,L)
    print(f"[{tag}] pontos {len(X)}: M1 (m_n>m_p) {int(m1.sum())}; M2 (m_L>m_p) {int(m2.sum())}; M1&M2 {int((m1&m2).sum())}; M1-M4 {int((m1&m2&m3&m4).sum())}; M1-M4&M5 {int((m1&m2&m3&m4&m5).sum())}")
    for nm_,mk in (("M1&M2",m1&m2),("M1-M4",m1&m2&m3&m4),("M1-M4&M5",m1&m2&m3&m4&m5)):
        if mk.sum(): 
            i=np.where(mk)[0][np.argmin(L[mk])]; print(f"   melhor L de raios entre {nm_}: {L[mk].min():.2f}; L<2: {int((L[mk]<2).sum())}; L<1: {int((L[mk]<1).sum())}; em (t,s,k,l)=({X[i][0]:.3f},{X[i][1]:.3f},{X[i][2]:.2f},{X[i][3]:.2f}) R=" + ",".join(f"{v:+.3f}" for v in R[i]))
        else: print(f"   {nm_}: nenhum ponto")
    print(f"   (controle) pontos com raios viaveis L<1 (sem massa): {int((L<1).sum())}; min L = {L.min():.2f}")
if __name__=="__main__":
    rng=np.random.default_rng(24); N=12000
    smp=qmc.LatinHypercube(d=4,seed=24).random(N)
    zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: TAB=np.array(P.map(work,list(X),chunksize=40))
    analyze(X,TAB,"fisico kappa>=0")
    X2=X.copy(); X2[:,2]=-X2[:,2]
    with Pool() as P: TAB2=np.array(P.map(work,list(X2),chunksize=40))
    analyze(X2,TAB2,"EXPLORATORIO nao fisico kappa<=0")
