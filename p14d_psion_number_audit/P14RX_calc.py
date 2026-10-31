import math, itertools, numpy as np
rng=np.random.default_rng(14)
Q9=np.array([1,-5,-5,10,4,4])/9
TH=np.array([math.pi*v/3 for v in range(6)]); XY=np.stack([np.cos(TH),np.sin(TH)],1)
PL=list(itertools.permutations(range(6),3))          # 120 colocacoes
TR={"p":(0,4,5),"n":(0,1,5),"Sm":(0,1,2),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
def arrays(t):
    q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum()
    A=np.zeros(120);B=np.zeros(120);C=np.zeros(120);R2=np.zeros(120)
    for i,st in enumerate(PL):
        A[i]=-sum(math.cos(TH[st[a]]-TH[t[a]]) for a in range(3))
        X=XY[list(st)]
        d=lambda a,b: np.linalg.norm(X[a]-X[b])
        B[i]=sum(q[a]*q[b]/d(a,b) for a in range(3) for b in range(a+1,3))
        C[i]=-sum(1/d(a,b) for a in range(3) for b in range(a+1,3) if ch[a]!=ch[b])
        c=X.mean(0); v=float(np.sum(q*np.sum((X-c)**2,1))); R2[i]=v/Qt if abs(Qt)>1e-12 else v
    return A,B,C,R2
AR={k:arrays(t) for k,t in TR.items()}
ks=np.arange(0,30.0001,0.25); GK,GL=np.meshgrid(ks,ks,indexing='ij'); GK=GK.ravel(); GL=GL.ravel(); G=len(GK)
def solve(key):
    A,B,C,R2=AR[key]; E=A[None,:]+GK[:,None]*B[None,:]+GL[:,None]*C[None,:]
    i=np.argmin(E,1); Es=np.sort(E,1); Emin=Es[:,0]
    amb=(Es[:,1]-Es[:,0]<1e-9)&(np.abs(R2[np.argsort(E,1)[:,1]]-R2[i])>1e-6)
    r=R2[i].copy(); r[amb]=np.nan; return Emin,r,i
S={k:solve(k) for k in TR}
Rrel=lambda k: S[k][1]/S["p"][1]
R={k:Rrel(k) for k in ("n","Sm","S0","L","Sp1","Sp2")}
# controle
i0=np.where((GK==0)&(GL==0))[0][0]
print(f"controle kappa=lambda=0: R_n = {R['n'][i0]:+.4f}, R_Sigma- = {R['Sm'][i0]:+.4f}")
def feasible(cn,cs0,cl,csp,w=(0.0164,0.15,0.125,0.175),use_sm=False):
    ok=(np.abs(R["n"]-cn)<w[0])&(np.abs(R["S0"]-cs0)<w[1])&(np.abs(R["L"]-cl)<w[2])&((np.abs(R["Sp1"]-csp)<w[3])|(np.abs(R["Sp2"]-csp)<w[3]))
    return ok
real=feasible(-0.164,0.0,0.125,0.975)
# alvo n a 10% relativo: janela 0,0164 ja usada
print(f"radii (n, Sigma0, Lambda, Sigma+): {int(real.sum())} de {G} pontos da grade viaveis")
okm=np.zeros(G,bool); info={}
Ep,En=S["p"][0],S["n"][0]; EL=S["L"][0]; ES0=S["S0"][0]
for cand in ("Sp1","Sp2"):
    ESp=S[cand][0]; dL=EL-Ep
    with np.errstate(all='ignore'):
        u=177.4/dL; sp=u*(ESp-Ep); s0=u*(ES0-EL)
        m=(dL>0)&((En-Ep)>0)&(np.abs(sp/251-1)<0.3)&(np.abs(s0/77-1)<0.3)
    info[cand]=(m,u,sp,s0)
    print(f"   massas com {cand}: {int((m&real).sum())} dos viaveis em raio passam sinal n-p e Sigma+-p, Sigma0-Lambda a 30%")
# lista de solucoes em raio
idx=np.where(real)[0]
if len(idx):
    print("   exemplos viaveis em raio (kappa, lambda): "+", ".join(f"({GK[i]:.2f},{GL[i]:.2f})" for i in idx[:8])+f" ... kappa em [{GK[idx].min():.2f},{GK[idx].max():.2f}], lambda em [{GL[idx].min():.2f},{GL[idx].max():.2f}]")
    i=idx[0]
    def plac(key):
        t=TR[key]; return PL[S[key][2][i]]
    print(f"   colocacoes no ponto ({GK[i]:.2f},{GL[i]:.2f}): "+"; ".join(f"{k}{TR[k]}->{plac(k)}" for k in ("p","n","S0","L","Sp1","Sp2")))
    print("   massas nesse ponto: "+", ".join(f"{c}: Sigma+-p={info[c][2][i]:+.0f} MeV, Sigma0-Lambda={info[c][3][i]:+.0f} MeV, n-p={(En[i]-Ep[i])*info[c][1][i] if np.isfinite(info[c][1][i]) else float('nan'):+.2f} MeV" for c in ("Sp1","Sp2")))
# teste de acaso
w=(0.0164,0.15,0.125,0.175); hits=0; N=2000
for _ in range(N):
    c=rng.uniform(-0.6,2.2,4)
    if feasible(c[0],c[1],c[2],c[3],w).any(): hits+=1
print(f"teste de acaso (alvos aleatorios, mesmas larguras): P_acaso = {hits/N:.3f} ({hits}/{N})")
# variante: apenas 3 alvos (sem n) e com n e Sigma-
for lbl,fn in (("com Sigma- (0,86+-0,21)",lambda: real&(np.abs(R['Sm']-0.86)<0.21)),):
    print(f"   {lbl}: {int(fn().sum())} pontos viaveis")
