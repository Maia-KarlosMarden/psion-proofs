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
if __name__=="__main__":
    rng=np.random.default_rng(26); N=16000
    smp=qmc.LatinHypercube(d=4,seed=26).random(N)
    zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,(smp[:,3]*2-1)*20],1)
    with Pool() as P: TAB=np.array(P.map(work,list(X),chunksize=40))
    np.savez('P14SGN_table.npz',X=X,TAB=TAB)
    R,T=unpack(TAB); L=Lrad(R); E={h:T[:,IDX[h],1] for h in H_}
    neg=X[:,3]<0; pos=X[:,3]>0
    print(f"pontos: {N}; lambda<0: {int(neg.sum())}; lambda>0: {int(pos.sum())}")
    # T1
    ok=mass_ok(E)&neg; D=E["L"]-E["n"]
    print(f"[T1 lambda<0] D>0: {int(((D>0)&neg).sum())}/{int(neg.sum())}; E_n>E_p: {int(((E['n']>E['p'])&neg).sum())}; D>0&R1&R2: {int(ok.sum())}")
    if ok.sum():
        i=np.where(ok)[0][np.argmin(L[ok])]; print(f"   melhor L de raios entre os que passam: {L[ok].min():.2f}; L<2: {int((L[ok]<2).sum())}; L<1: {int((L[ok]<1).sum())}; em (t,s,k,l)=({X[i][0]:.3f},{X[i][1]:.3f},{X[i][2]:.2f},{X[i][3]:.2f}) R="+",".join(f"{v:+.3f}" for v in R[i]))
    print(f"   (controle) lambda<0: min L de raios (sem massa) = {L[neg].min():.2f}; L<1: {int(((L<1)&neg).sum())}")
    # T2
    ws=np.array([-3,-2,-1,-0.5,0,0.5,1,2,3.]); W=np.array(np.meshgrid(ws,ws,ws,indexing='ij')).reshape(3,-1).T   # (729,3)
    def anyw(r1t,r2t,idxpts):
        out=np.zeros(len(idxpts),bool); wbest=[None]*len(idxpts)
        for n_,i in enumerate(idxpts):
            v={h:W@T[i,IDX[h],2:5] for h in H_}   # (729,) massas por peso
            m=mass_ok(v,r1t,r2t)
            if m.any(): out[n_]=True
        return out
    ip=np.where(pos)[0]; a_real=anyw(0.437,1.425,ip)
    print(f"[T2 lambda>0] pontos com ALGUM peso passando (alvo real): {int(a_real.sum())}/{len(ip)} ({100*a_real.mean():.1f}%)")
    sub=rng.choice(len(ip),1500,replace=False)
    fr=[]
    for _ in range(40):
        r1t=rng.uniform(0.1,1.0); r2t=rng.uniform(0.5,3.0); fr.append(anyw(r1t,r2t,ip[sub]).mean())
    fr=np.array(fr); real_sub=a_real[sub].mean()
    print(f"   controle de vacuidade: alvos de massa aleatorios (40, 1500 pontos): fracao media {100*fr.mean():.1f}% [min {100*fr.min():.1f}, max {100*fr.max():.1f}]; alvo real nessa amostra {100*real_sub:.1f}%")
    if a_real.sum():
        idx=ip[a_real]; Lm=L[idx]; j=np.argmin(Lm)
        print(f"   entre os que passam em massa: melhor L de raios = {Lm.min():.2f}; L<2: {int((Lm<2).sum())}; L<1: {int((Lm<1).sum())}; em (t,s,k,l)=({X[idx[j]][0]:.3f},{X[idx[j]][1]:.3f},{X[idx[j]][2]:.2f},{X[idx[j]][3]:.2f})")
        print(f"   (controle) lambda>0 sem massa: min L = {L[pos].min():.2f}; L<1: {int(((L<1)&pos).sum())}")
