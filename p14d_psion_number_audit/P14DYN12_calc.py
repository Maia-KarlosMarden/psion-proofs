import math, itertools, numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
rng=np.random.default_rng(19)
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
tv=[0,0.1,0.3,1,3]; sv=[0,0.1,0.3,1,3]; kv=[0,2,8]; lv=[0,2,8]
grid=[(t,s,k,l) for t in tv for s in sv for k in kv for l in lv]
Rtab={h:np.zeros(len(grid)) for h in TR}
def gs(h,t,s,k,l):
    A,B,C,R2=AR[h]; dg=A+k*B+l*C
    if t==0 and s==0:
        o=np.sort(dg)[:2]
        if o[1]-o[0]<1e-9:
            i0=np.where(dg<o[0]+1e-9)[0]
            return float(R2[i0[0]]) if np.ptp(R2[i0])<1e-6 else np.nan
        return float(R2[np.argmin(dg)])
    H=sp.diags(dg)-t*T1-s*S1
    w,v=eigsh(H.tocsc(),k=2,which='SA')
    o=np.argsort(w); w=w[o]; v=v[:,o]
    if w[1]-w[0]<1e-8:
        r0=(v[:,0]**2*R2).sum(); r1=(v[:,1]**2*R2).sum()
        if abs(r0-r1)>1e-6: return np.nan
    return float((v[:,0]**2*R2).sum())
for g,(t,s,k,l) in enumerate(grid):
    for h in TR: Rtab[h][g]=gs(h,t,s,k,l)
rp=Rtab["p"]; R={h:Rtab[h]/rp for h in TR if h!="p"}
W=(0.0164,0.15,0.10,0.175)
def feas(cn,cs0,cl,csp):
    return (np.abs(R["n"]-cn)<W[0])&(np.abs(R["S0"]-cs0)<W[1])&(np.abs(R["L"]-cl)<W[2])&((np.abs(R["Sp1"]-csp)<W[3])|(np.abs(R["Sp2"]-csp)<W[3]))
real=feas(-0.164,0.0,0.15,0.975); G=len(grid)
print(f"raios completos (n, Sigma0, Lambda, Sigma+): {int(real.sum())} de {G}")
for nm,m in (("n",np.abs(R['n']+0.164)<W[0]),("Sigma0",np.abs(R['S0'])<W[1]),("Lambda",np.abs(R['L']-0.15)<W[2]),("Sigma+",(np.abs(R['Sp1']-0.975)<W[3])|(np.abs(R['Sp2']-0.975)<W[3]))):
    print(f"   so {nm}: {int(m.sum())} pontos")
for i in np.where(real)[0][:5]:
    print("   viavel (t,s,kappa,lambda)=",grid[i],{h:round(float(R[h][i]),3) for h in R})
L=R['L'][~np.isnan(R['L'])]; print(f"   R_L: min {L.min():+.3f} max {L.max():+.3f}; R_n: min {np.nanmin(R['n']):+.3f} max {np.nanmax(R['n']):+.3f}")
hits=0;N=2000
for _ in range(N):
    c=rng.uniform(-0.6,2.2,4)
    if feas(*c).any(): hits+=1
print(f"teste de acaso: P_acaso = {hits/N:.3f} ({hits}/{N})")
