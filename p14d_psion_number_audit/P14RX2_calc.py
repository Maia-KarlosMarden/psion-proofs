import math, itertools, numpy as np
rng=np.random.default_rng(15)
Q9=np.array([1,-5,-5,10,4,4])/9; NS=12
TH=np.array([2*math.pi*s/NS for s in range(NS)]); XY=np.stack([np.cos(TH),np.sin(TH)],1)
PL=list(itertools.permutations(range(NS),3)); PA=np.array(PL)
TR={"p":(0,4,5),"n":(0,1,5),"Sm":(0,1,2),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def arrays(t):
    q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum()
    X=XY[PA]                                        # (N,3,2)
    A=-sum(np.cos(TH[PA[:,a]]-TH[2*t[a]]) for a in range(3))
    d=lambda a,b: np.linalg.norm(X[:,a]-X[:,b],axis=1)
    B=sum(q[a]*q[b]/d(a,b) for a in range(3) for b in range(a+1,3))
    Cc=-sum(1/d(a,b) for a in range(3) for b in range(a+1,3) if ch[a]!=ch[b]) if any(ch[a]!=ch[b] for a in range(3) for b in range(a+1,3)) else np.zeros(len(PL))
    c=X.mean(1,keepdims=True); v=(q[None,:]*((X-c)**2).sum(2)).sum(1); R2=v/Qt if abs(Qt)>1e-12 else v
    R2=np.where(np.abs(R2)<1e-12,0.0,R2)
    return A,B,Cc,R2
AR={k:arrays(t) for k,t in TR.items()}
ks=np.arange(0,30.0001,0.5); GK,GL=np.meshgrid(ks,ks,indexing='ij'); GK=GK.ravel(); GL=GL.ravel(); G=len(GK)
def solve(key):
    A,B,C,R2=AR[key]; E=A[None,:]+GK[:,None]*B[None,:]+GL[:,None]*C[None,:]
    o=np.argsort(E,1)[:,:2]; i=o[:,0]; e0=E[np.arange(G),i]; e1=E[np.arange(G),o[:,1]]
    amb=(e1-e0<1e-9)&(np.abs(R2[o[:,1]]-R2[i])>1e-6); r=R2[i].copy(); r[amb]=np.nan; return e0,r,i
S={k:solve(k) for k in TR}
R={k:S[k][1]/S["p"][1] for k in ("n","Sm","S0","L","Sp1","Sp2")}
i0=int(np.where((GK==0)&(GL==0))[0][0]); print(f"controle kappa=lambda=0: R_n = {R['n'][i0]:+.4f}, R_Sigma- = {R['Sm'][i0]:+.4f}")
W=(0.0164,0.15,0.10,0.175)
def feas(cn,cs0,cl,csp,w=W):
    return (np.abs(R["n"]-cn)<w[0])&(np.abs(R["S0"]-cs0)<w[1])&(np.abs(R["L"]-cl)<w[2])&((np.abs(R["Sp1"]-csp)<w[3])|(np.abs(R["Sp2"]-csp)<w[3]))
real=feas(-0.164,0.0,0.15,0.975)       # Lambda: centro 0,15, meia-largura 0,10 -> [0,05; 0,25]
print(f"fases distintas (colocacoes) na grade: {len(set(tuple(S[k][2][i] for k in TR) for i in range(G)))}")
print(f"raios (n, Sigma0, Lambda, Sigma+): {int(real.sum())} de {G} pontos viaveis")
Ep,En,EL,ES0=S["p"][0],S["n"][0],S["L"][0],S["S0"][0]; ok_m={}
for cand in ("Sp1","Sp2"):
    ESp=S[cand][0]; dL=EL-Ep
    with np.errstate(all='ignore'):
        u=177.4/dL; sp=u*(ESp-Ep); s0=u*(ES0-EL)
        m=(dL>0)&((En-Ep)>0)&(np.abs(sp/251-1)<0.3)&(np.abs(s0/77-1)<0.3)
    ok_m[cand]=(m,u,sp,s0)
    print(f"   massas ({cand}): {int((m&real).sum())} dos viaveis em raio passam")
idx=np.where(real)[0]
if len(idx):
    print(f"   kappa em [{GK[idx].min():.1f},{GK[idx].max():.1f}], lambda em [{GL[idx].min():.1f},{GL[idx].max():.1f}]")
    for i in idx[:3]:
        print(f"   ({GK[i]:.1f},{GL[i]:.1f}): "+"; ".join(f"{k}{TR[k]}->{tuple(PL[S[k][2][i]])}" for k in ("p","n","S0","L","Sp1","Sp2")))
        print("      R: "+", ".join(f"{k}={R[k][i]:+.3f}" for k in ("n","Sm","S0","L","Sp1","Sp2")))
        for c in ("Sp1","Sp2"): print(f"      massas ({c}): Sigma+-p = {ok_m[c][2][i]:+.0f} MeV, Sigma0-Lambda = {ok_m[c][3][i]:+.0f} MeV, n-p>0: {bool(En[i]-Ep[i]>0)}")
hits=0;N=2000
for _ in range(N):
    c=rng.uniform(-0.6,2.2,4)
    if feas(c[0],c[1],c[2],c[3]).any(): hits+=1
print(f"teste de acaso: P_acaso = {hits/N:.3f} ({hits}/{N})")
print(f"   com Sigma- (0,86+-0,21) adicionalmente: {int((real&(np.abs(R['Sm']-0.86)<0.21)).sum())} pontos")
