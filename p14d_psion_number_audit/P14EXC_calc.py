import math, itertools, numpy as np
rng=np.random.default_rng(17)
Q9=np.array([1,-5,-5,10,4,4])/9
TH=np.array([math.pi/3*s for s in range(6)]); XY=np.stack([np.cos(TH),np.sin(TH)],1)
TR={"p":(0,4,5),"n":(0,1,5),"Sm":(0,1,2),"S0":(0,2,4),"L":(1,2,3),"Sp1":(1,3,5),"Sp2":(2,3,4)}
ks=np.arange(0,30.0001,0.5); GK,GL=np.meshgrid(ks,ks,indexing='ij'); GK=GK.ravel(); GL=GL.ravel(); G=len(GK)
def perms(t):
    out=[]
    for P in itertools.permutations(t):      # P[a] = vertice ocupado pelo rotulo t[a]
        X=XY[list(P)]; q=Q9[list(t)]; ch=['L' if k<=2 else 'R' for k in t]; Qt=q.sum()
        d=lambda a,b: np.linalg.norm(X[a]-X[b])
        A=-sum(math.cos(TH[P[a]]-TH[t[a]]) for a in range(3))
        B=sum(q[a]*q[b]/d(a,b) for a in range(3) for b in range(a+1,3))
        C=-sum(1/d(a,b) for a in range(3) for b in range(a+1,3) if ch[a]!=ch[b])
        c=X.mean(0); v=sum(q[a]*((X[a]-c)**2).sum() for a in range(3)); R2=v/Qt if abs(Qt)>1e-12 else v
        out.append((A,B,C,0.0 if abs(R2)<1e-12 else R2,P))
    return out
PM={k:perms(t) for k,t in TR.items()}
def solve(k):
    A=np.array([x[0] for x in PM[k]]);B=np.array([x[1] for x in PM[k]]);C=np.array([x[2] for x in PM[k]]);R2=np.array([x[3] for x in PM[k]])
    E=A[None,:]+GK[:,None]*B[None,:]+GL[:,None]*C[None,:]; o=np.argsort(E,1)[:,:2]; i=o[:,0]
    e0=E[np.arange(G),i]; e1=E[np.arange(G),o[:,1]]; r=R2[i].copy(); amb=(e1-e0<1e-9)&(np.abs(R2[o[:,1]]-R2[i])>1e-6); r[amb]=np.nan
    return e0,r,i
S={k:solve(k) for k in TR}
print("valores de R2 acessiveis por hadron (todas as 6 permutacoes):")
rp_all=sorted(set(round(x[3],4) for x in PM["p"])); print("  p:",rp_all)
for k in ("n","S0","L","Sp1","Sp2"): print(f"  {k}:",sorted(set(round(x[3],4) for x in PM[k])))
R={k:S[k][1]/S["p"][1] for k in ("n","Sm","S0","L","Sp1","Sp2")}
W=(0.0164,0.15,0.10,0.175)
def feas(cn,cs0,cl,csp):
    return (np.abs(R["n"]-cn)<W[0])&(np.abs(R["S0"]-cs0)<W[1])&(np.abs(R["L"]-cl)<W[2])&((np.abs(R["Sp1"]-csp)<W[3])|(np.abs(R["Sp2"]-csp)<W[3]))
real=feas(-0.164,0.0,0.15,0.975)
print(f"raios completos (n, Sigma0, Lambda, Sigma+): {int(real.sum())} de {G}")
for nm,m in (("n",np.abs(R['n']+0.164)<W[0]),("Sigma0",np.abs(R['S0'])<W[1]),("Lambda",np.abs(R['L']-0.15)<W[2]),("Sigma+",(np.abs(R['Sp1']-0.975)<W[3])|(np.abs(R['Sp2']-0.975)<W[3]))):
    print(f"   so {nm}: {int(m.sum())} pontos")
print("   valores de R_L alcancados:",sorted(set(np.round(R['L'][~np.isnan(R['L'])],3))))
print("   valores de R_n alcancados:",sorted(set(np.round(R['n'][~np.isnan(R['n'])],3))))
hits=0;N=2000
for _ in range(N):
    c=rng.uniform(-0.6,2.2,4)
    if feas(*c).any(): hits+=1
print(f"teste de acaso: P_acaso = {hits/N:.3f} ({hits}/{N})")
