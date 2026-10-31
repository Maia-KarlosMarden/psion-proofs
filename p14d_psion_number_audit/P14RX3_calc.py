import io, contextlib, runpy, math, numpy as np
rng=np.random.default_rng(16)
with contextlib.redirect_stdout(io.StringIO()): X=runpy.run_path('P14RX2_calc.py',run_name='lib')
AR,PL,PA,TR,XY,TH,S,GK,GL,G,Q9=X['AR'],X['PL'],X['PA'],X['TR'],X['XY'],X['TH'],X['S'],X['GK'],X['GL'],X['G'],X['Q9']
bs=np.arange(0,30.0001,0.5); NB=len(bs)
def Dvec(t):
    if 2 in t and 3 in t:
        a,b=t.index(2),t.index(3); return -1/np.linalg.norm(XY[PA[:,a]]-XY[PA[:,b]],axis=1)
    return np.zeros(len(PL))
D={k:Dvec(t) for k,t in TR.items()}
def solveb(key):
    A,B,C,R2=AR[key]; E0=A[None,:]+GK[:,None]*B[None,:]+GL[:,None]*C[None,:]
    e=np.zeros((NB,G)); r=np.zeros((NB,G))
    for j,b in enumerate(bs):
        E=E0+b*D[key][None,:]; ip=np.argpartition(E,1,axis=1)[:,:2]
        v=np.take_along_axis(E,ip,1); o=np.argsort(v,1); i1=np.take_along_axis(ip,o[:,:1],1)[:,0]; i2=np.take_along_axis(ip,o[:,1:2],1)[:,0]
        e1=E[np.arange(G),i1]; e2=E[np.arange(G),i2]; rr=R2[i1].copy(); amb=(e2-e1<1e-9)&(np.abs(R2[i2]-R2[i1])>1e-6); rr[amb]=np.nan
        e[j]=e1; r[j]=rr
    return e,r
EL,RL2=solveb("L"); ES2,RS2=solveb("Sp2")
rp=S["p"][1]; Rn=S["n"][1]/rp; R0=S["S0"][1]/rp; Rs1=S["Sp1"][1]/rp
RL=RL2/rp[None,:]; Rs2=RS2/rp[None,:]
W=(0.0164,0.15,0.10,0.175)
def feas(cn,cs0,cl,csp):
    return (np.abs(Rn-cn)<W[0])[None,:]&(np.abs(R0-cs0)<W[1])[None,:]&(np.abs(RL-cl)<W[2])&((np.abs(Rs1-csp)<W[3])[None,:]|(np.abs(Rs2-csp)<W[3]))
real=feas(-0.164,0.0,0.15,0.975)
print(f"raios (n, Sigma0, Lambda, Sigma+): {int(real.sum())} de {NB*G} pontos (kappa, lambda, beta) viaveis; beta em "
      +(f"[{bs[np.where(real.any(1))[0]].min():.1f},{bs[np.where(real.any(1))[0]].max():.1f}]" if real.any() else "-"))
Ep,En,E0s=S["p"][0],S["n"][0],S["S0"][0]; ESp1=S["Sp1"][0]
dL=EL-Ep[None,:]
with np.errstate(all='ignore'):
    u=177.4/dL
    ok={}
    for cand,ESp in (("Sp1",np.broadcast_to(ESp1,(NB,G))),("Sp2",ES2)):
        sp=u*(ESp-Ep[None,:]); s0=u*(E0s[None,:]-EL)
        ok[cand]=(dL>0)&((En-Ep)[None,:]>0)&(np.abs(sp/251-1)<0.3)&(np.abs(s0/77-1)<0.3)
    mass_pass=(ok["Sp1"]|ok["Sp2"])&real
    # diagnostico: onde falha
    sign_np=((En-Ep)>0)[None,:]&real; dLpos=(dL>0)&real
print(f"   entre os viaveis em raio: sinal m_n > m_p em {int(sign_np.sum())}; m_Lambda > m_p em {int(dLpos.sum())}; ambas as massas a 30% em {int(mass_pass.sum())}")
idx=np.argwhere(real)
if len(idx):
    for j,i in idx[::max(1,len(idx)//4)][:4]:
        print(f"   (kappa, lambda, beta) = ({GK[i]:.1f}, {GL[i]:.1f}, {bs[j]:.1f}): R_n={Rn[i]:+.3f} R_S0={R0[i]:+.3f} R_L={RL[j,i]:+.3f} R_Sp1={Rs1[i]:+.3f} R_Sp2={Rs2[j,i]:+.3f} | E_n-E_p={En[i]-Ep[i]:+.3f}, E_L-E_p={dL[j,i]:+.3f}, E_Sp2-E_p={ES2[j,i]-Ep[i]:+.3f}")
hits=0;N=2000
for _ in range(N):
    c=rng.uniform(-0.6,2.2,4)
    if feas(*c).any(): hits+=1
print(f"teste de acaso: P_acaso = {hits/N:.3f} ({hits}/{N})")
print(f"   com Sigma- (0,86+-0,21): {int((real&(np.abs(S['Sm'][1]/rp-0.86)<0.21)[None,:]).sum())} pontos")
