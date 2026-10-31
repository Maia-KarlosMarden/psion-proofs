import math, itertools, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
TH=[math.pi*v/3 for v in range(6)]
XY=np.array([[math.cos(a),math.sin(a)] for a in TH])
def build(t_labels):
    S=[p for p in itertools.permutations(range(6),3)]            # S[i] = (v1,v2,v3): vertice de cada rotulo
    idx={s:i for i,s in enumerate(S)}
    return S,idx
def hamiltonian(tl,eps,t,s):
    S,idx=build(tl); n=len(S); H=np.zeros((n,n))
    for i,st in enumerate(S):
        H[i,i]=-eps*sum(math.cos(TH[st[a]]-TH[tl[a]]) for a in range(3))
        occ=set(st)
        for a in range(3):
            for d,amp in ((1,t),(-1,t),(3,s)):
                if amp==0: continue
                w=(st[a]+d)%6
                if w in occ: continue
                new=list(st); new[a]=w; H[idx[tuple(new)],i]-=amp
        if s!=0:                                                  # troca de posicoes de dois psions antipodais
            for a in range(3):
                for b in range(a+1,3):
                    if (st[a]-st[b])%6==3:
                        new=list(st); new[a],new[b]=st[b],st[a]; H[idx[tuple(new)],i]-=s
    return H,S
def r2state(tl,st):
    q=Q9[list(tl)]; X=XY[list(st)]; c=X.mean(0); return float(np.sum(q*np.sum((X-c)**2,1)))
def gs(tl,eps,t,s):
    H,S=hamiltonian(tl,eps,t,s); H=(H+H.T)/2
    w,v=np.linalg.eigh(H); g=v[:,0]; P=g**2
    nat=S.index(tuple(tl))
    r2=sum(P[i]*r2state(tl,S[i]) for i in range(len(S))); Q=Q9[list(tl)].sum()
    return dict(E0=w[0],gap=w[1]-w[0],w=w,pnat=P[nat],r2=r2/Q if abs(Q)>1e-12 else r2)
T={"p":(0,4,5),"n":(0,1,5),"Sigma-":(0,1,2),"Sigma0":(0,2,4),"Lambda":(1,2,3),"Sigma+ (2,3,4)":(2,3,4)}
eps=1.0
# controle
c=gs(T["n"],1,0,0); cp=gs(T["p"],1,0,0); cs=gs(T["Sigma-"],1,0,0)
print(f"controle t=s=0: E0(n) = {c['E0']:.4f} (esperado -3), R_n = {c['r2']/cp['r2']:+.4f}, R_Sigma- = {cs['r2']/cp['r2']:+.4f}")
print(f"limite delocalizado (eps=0, t=1): R_n = {gs(T['n'],0,1,0)['r2']/gs(T['p'],0,1,0)['r2']:+.4f}, R_Sigma- = {gs(T['Sigma-'],0,1,0)['r2']/gs(T['p'],0,1,0)['r2']:+.4f}")
print("custos de excitacao de um psion (eps): passo 60 graus =",round(1-math.cos(math.pi/3),3),"; 120 graus =",round(1-math.cos(2*math.pi/3),3),"; antipoda =",2.0)
for sname,sfac in (("s=0",0.0),("s=t",1.0)):
    print(f"--- {sname} ---")
    print("  t/eps   E0(n)   dE0/(t^2/eps)  gap(n)  p_nat(n)  p_nat(p) |  R_n    R_Sig-  R_Sig0  R_Lam   R_Sig+(234)")
    for t in (0,0.02,0.05,0.1,0.2,0.5,1,2,5,20):
        s=sfac*t; rp=gs(T["p"],eps,t,s)
        row={k:gs(v,eps,t,s) for k,v in T.items()}
        dE=(row["n"]["E0"]+3)/(t*t) if t>0 else float('nan')
        print(f"  {t:5.2f}  {row['n']['E0']:+7.3f}  {dE:+10.2f}   {row['n']['gap']:6.3f}  {row['n']['pnat']:7.3f}  {rp['pnat']:7.3f} | "
              +"  ".join(f"{row[k]['r2']/rp['r2']:+6.3f}" for k in ("n","Sigma-","Sigma0","Lambda","Sigma+ (2,3,4)")))
# espectro de baixa energia (n, t/eps = 0.2, s=0 e s=t)
for sfac in (0,1):
    H,S=hamiltonian(T["n"],1,0.2,0.2*sfac); w=np.linalg.eigvalsh((H+H.T)/2)
    print(f"espectro n, t/eps=0,2, s={'t' if sfac else '0'}: 8 niveis mais baixos = {np.round(w[:8],3).tolist()}; largura total = {w[-1]-w[0]:.2f}")
