import io, contextlib, math, runpy, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    N_=runpy.run_path('P14SIM_network.py',run_name='lib')
one,OJS,J,PL,PR,V3=N_['one'],N_['OJS'],N_['J'],N_['PL'],N_['PR'],N_['V3']
mk=np.array([k-2.5 for k in range(6)]); Jp=J[0]+1j*J[1]
X=one(Jp@PL)@one(Jp@PR); H0d=one(np.diag(np.abs(mk)).astype(complex))
bond=np.zeros((6,6),complex); bond[2,3]=bond[3,2]=0.5; Bond=one(bond)
proj=lambda M: V3.conj()@M@V3.T
H0=proj(H0d+OJS); Xs=proj(X); Bs=proj(Bond)
th=2*math.pi*np.arange(96)/96
def Ecurve(extra):
    return np.array([np.linalg.eigvalsh(H0+(1/6)*(np.exp(1j*t)*Xs+np.exp(-1j*t)*Xs.conj().T)+extra)[0] for t in th])
def harm(E): return np.abs(np.fft.rfft(E))/96*2
dC=np.array([0,15,24,27,24,15],float); s=np.arange(6)
print("T2: DeltaC(s) em cos(n pi s/3): n=0..3:",[round(float(np.sum(dC*np.cos(n*np.pi*s/3))*(1 if n in(0,3) else 2)/6),3) for n in range(4)])
tgt=dC/27
cands={"Bond d=0.05":0.05*Bs}
eta=0.02
for j in range(1,7):
    Xj=np.linalg.matrix_power(Xs,j); cands[f"X^{j} eta=0.02"]=eta*(Xj+Xj.conj().T); print(f"||X^{j}|| no setor = {np.abs(Xj).max():.3e}")
for nm,ex in cands.items():
    if np.abs(ex).max()<1e-12: print(f"{nm}: operador nulo no setor de 3 psions"); continue
    E=Ecurve(ex); F=harm(E); idx=[int(round(k*96/6)) for k in range(6)]; v=E[idx]-E[idx[0]]
    nv=v/np.abs(v).max() if np.abs(v).max()>1e-14 else v
    dev=np.abs(nv-tgt).max() if np.abs(v).max()>1e-14 else np.nan
    print(f"{nm}: harmonicas n=1..6 = {[f'{F[n]:.1e}' for n in range(1,7)]}; forma E(theta_k)-E(0) normalizada = {np.round(nv,3)}; desvio max vs DeltaC/27 = {dev:.3f}")
