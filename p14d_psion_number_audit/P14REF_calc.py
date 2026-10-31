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

def evalz(z):
    t,s,k,l,a,b,g=z
    r={h:gse(h,t,s,k,l) for h in H_}
    if any(np.isnan(r[h][0]) for h in H_): return 50.0,None
    rp=r["p"][0]; R=[r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]
    Lr=max(abs(R[0]+0.164)/0.0164,abs(R[1])/0.15,abs(R[2]-0.15)/0.10,min(abs(R[3]-0.975),abs(R[4]-0.975))/0.175)
    M={h:a*r[h][2]+b*r[h][3]+g*r[h][4] for h in H_}
    D=M["L"]-M["n"]
    if D<=0: Lm=5+abs(D)
    else:
        v1=abs((M["S0"]-M["L"])/D/0.437-1)/0.3; v2=min(abs((M["Sp1"]-M["p"])/D/1.425-1),abs((M["Sp2"]-M["p"])/D/1.425-1))/0.3; Lm=max(v1,v2)
    return max(Lr,Lm),(R,D,(M["S0"]-M["L"])/D if D>0 else np.nan,Lr,Lm)
def job(arg):
    var,z0=arg
    lo=np.array([0,0,0,0]+([-3]*3 if var=="A" else [0]*3),float); hi=np.array([6,6,20,20]+([3]*3 if var=="A" else [1.5]*3),float)
    f=lambda z: evalz(np.clip(z,lo,hi))[0]
    r=minimize(f,np.clip(z0,lo,hi),method='Nelder-Mead',options={'maxiter':250,'xatol':1e-4,'fatol':1e-5})
    z=np.clip(r.x,lo,hi); return var,float(r.fun),z
if __name__=="__main__":
    d=np.load('P14SGN_table.npz'); X=d['X']; TAB=d['TAB']
    R,T=unpack(TAB); L=Lrad(R); pos=np.where(X[:,3]>0)[0]
    ws=np.array([-3,-2,-1,-0.5,0,0.5,1,2,3.]); W=np.array(np.meshgrid(ws,ws,ws,indexing='ij')).reshape(3,-1).T
    wsB=np.array([0,0.5,1,1.5]); WB=np.array(np.meshgrid(wsB,wsB,wsB,indexing='ij')).reshape(3,-1).T
    def best_start(i,Ws):
        v={h:Ws@T[i,IDX[h],2:5] for h in H_}; D=v["L"]-v["n"]
        with np.errstate(all='ignore'):
            v1=np.abs((v["S0"]-v["L"])/D/0.437-1)/0.3; v2=np.minimum(np.abs((v["Sp1"]-v["p"])/D/1.425-1),np.abs((v["Sp2"]-v["p"])/D/1.425-1))/0.3
        Lm=np.where(D<=0,5+np.abs(D),np.maximum(v1,v2)); Lc=np.maximum(L[i],Lm); j=np.argmin(Lc); return Lc[j],Ws[j]
    starts=[]
    for var,Ws in (("A",W),("B",WB)):
        sc=[]
        for i in pos:
            lc,w=best_start(i,Ws); sc.append((lc,i,w))
        sc.sort(key=lambda a:a[0])
        print(f"variante {var}: melhor L_c na tabela = {sc[0][0]:.2f}")
        for lc,i,w in sc[:12]: starts.append((var,np.concatenate([X[i],w])))
        starts.append((var,np.array([0.001,0.069,0.26,0.11,1,1,1.])))
    with Pool() as P: res=P.map(job,starts)
    for var in ("A","B"):
        rr=sorted([r for r in res if r[0]==var],key=lambda a:a[1]); print(f"variante {var}: min L_c apos refino = {rr[0][1]:.3f} ({'VIAVEL' if rr[0][1]<1 else 'inviavel'}); viaveis: {sum(1 for r in rr if r[1]<1)}/{len(rr)}")
        for _,f,z in rr[:3]:
            Lc,inf=evalz(z); Rr,D,r1,Lr,Lm=inf
            nat=bool(np.all(np.abs(z[4:])<=1.5) and (np.all(z[4:]>=0) or np.all(z[4:]<=0)))
            print(f"   L_c={f:.3f} (L_r={Lr:.2f}, L_m={Lm:.2f}) (t,s,k,l)=({z[0]:.3f},{z[1]:.3f},{z[2]:.2f},{z[3]:.2f}) pesos (a,b,g)=({z[4]:+.2f},{z[5]:+.2f},{z[6]:+.2f}) natural={nat}; R="+",".join(f"{v:+.3f}" for v in Rr)+f"; D={D:.3f}")
