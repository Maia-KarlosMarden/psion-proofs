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
tv=[0,0.03,0.06,0.1,0.15,0.2,0.3,0.5,0.8,1.2,2,3.5]; sv=tv; kv=[0,0.5,1,2,3,5,8,12]; lv=kv
grid=[(t,s,k,l) for t in tv for s in sv for k in kv for l in lv]
def gs(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C
    if t==0 and s==0:
        o=np.sort(dg)[:2]
        if o[1]-o[0]<1e-9:
            i0=np.where(dg<o[0]+1e-9)[0]
            return (float(R2[i0[0]]),float(o[0])) if np.ptp(R2[i0])<1e-6 else (np.nan,np.nan)
        return float(R2[np.argmin(dg)]),float(o[0])
    H=sp.diags(dg)-t*T1-s*S1
    w,v=eigsh(H.tocsc(),k=2,which='SA'); o=np.argsort(w); w=w[o]; v=v[:,o]
    if w[1]-w[0]<1e-8 and abs((v[:,0]**2*R2).sum()-(v[:,1]**2*R2).sum())>1e-6: return np.nan,np.nan
    return float((v[:,0]**2*R2).sum()),float(w[0])
def work(g):
    t,s,k,l=g; return [gs(h,t,s,k,l) for h in TR]
if __name__=="__main__":
    with Pool() as P: out=P.map(work,grid,chunksize=40)
    out=np.array(out); H=list(TR); G=len(grid)
    R2t={h:out[:,i,0] for i,h in enumerate(H)}; Et={h:out[:,i,1] for i,h in enumerate(H)}
    R={h:R2t[h]/R2t["p"] for h in H if h!="p"}
    W=(0.0164,0.15,0.10,0.175)
    def feas(cn,cs0,cl,csp):
        return (np.abs(R["n"]-cn)<W[0])&(np.abs(R["S0"]-cs0)<W[1])&(np.abs(R["L"]-cl)<W[2])&((np.abs(R["Sp1"]-csp)<W[3])|(np.abs(R["Sp2"]-csp)<W[3]))
    real=feas(-0.164,0.0,0.15,0.975)
    print(f"raios completos (n, Sigma0, Lambda, Sigma+): {int(real.sum())} de {G}")
    for nm,m in (("n",np.abs(R['n']+0.164)<W[0]),("Sigma0",np.abs(R['S0'])<W[1]),("Lambda",np.abs(R['L']-0.15)<W[2]),("Sigma+",(np.abs(R['Sp1']-0.975)<W[3])|(np.abs(R['Sp2']-0.975)<W[3]))):
        print(f"   so {nm}: {int(m.sum())} pontos")
    nl=(np.abs(R['n']+0.164)<W[0])&(np.abs(R['L']-0.15)<W[2]); print(f"   n e Lambda juntos: {int(nl.sum())}")
    nls=nl&((np.abs(R['Sp1']-0.975)<W[3])|(np.abs(R['Sp2']-0.975)<W[3])); print(f"   n, Lambda e Sigma+ juntos: {int(nls.sum())}")
    dL=Et["L"]-Et["p"]
    with np.errstate(all='ignore'):
        u=177.4/dL; res={}
        for c in ("Sp1","Sp2"):
            sp_=u*(Et[c]-Et["p"]); s0=u*(Et["S0"]-Et["L"])
            res[c]=(dL>0)&((Et["n"]-Et["p"])>0)&(np.abs(sp_/251-1)<0.3)&(np.abs(s0/77-1)<0.3)
    ms=(res["Sp1"]|res["Sp2"])&real
    print(f"   entre os viaveis em raio: massas (n>p, L>p, Sigma+ e Sigma0-Lambda a 30%): {int(ms.sum())}")
    for i in np.where(real)[0][:6]:
        print("   viavel (t,s,kappa,lambda)=",grid[i],{h:round(float(R[h][i]),3) for h in R},"dE_n=%.3f dE_L=%.3f"%(Et['n'][i]-Et['p'][i],dL[i]))
    rng=np.random.default_rng(20); hits=0;N=2000
    for _ in range(N):
        c=rng.uniform(-0.6,2.2,4)
        if feas(*c).any(): hits+=1
    print(f"teste de acaso (grade densa): P_acaso = {hits/N:.3f} ({hits}/{N})")
