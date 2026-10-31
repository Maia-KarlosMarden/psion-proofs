# ADENDO POS-HOC (nao pre-registrado antes de ver o ponto nominal viavel de V0 na tabela de HYB2): localizar o ponto e checar detalhes/estabilidade.
import io, contextlib, runpy, numpy as np
from multiprocessing import Pool
from scipy.stats import qmc
with contextlib.redirect_stdout(io.StringIO()): H=runpy.run_path('P14HYB2_calc.py',run_name='lib')
full,Lc,pert,A=H['full'],H['Lc'],H['pert'],H['A']; DIRS,MX,H_=A['DIRS'],A['MX'],A['H_']
def f(x): return Lc(x,"V0")
if __name__=="__main__":
    rng=np.random.default_rng(43); N=8000; smp=qmc.LatinHypercube(d=4,seed=43).random(N); zt=rng.random(N)<0.05; zs=rng.random(N)<0.05
    lg=lambda u: np.exp(np.log(0.01)+u*(np.log(6)-np.log(0.01)))
    X=np.stack([np.where(zt,0,lg(smp[:,0])),np.where(zs,0,lg(smp[:,1])),smp[:,2]*20,smp[:,3]*20],1)
    with Pool() as P: L=np.array(P.map(f,list(X),chunksize=80))
    idx=np.argsort(L)[:3]; print("3 menores L_c (V0) na tabela:",[(round(float(L[i]),3),np.round(X[i],3).tolist()) for i in idx])
    z=X[idx[0]]; Lr,Av,R=full(z)
    print(f"ponto {np.round(z,4)}: L_r={Lr:.3f}; R_n={R[0]:+.3f} R_S0={R[1]:+.3f} R_L={R[2]:+.3f} R_Sp1={R[3]:+.3f} R_Sp2={R[4]:+.3f} R_Sm={R[5]:+.3f}")
    M={h:DIRS@Av[h] for h in Av}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        r1=(M["S0"]-M["L"])/D; a=(M["Sp1"]-M["p"])/D; b=(M["Sp2"]-M["p"])/D; q=D/MX
        ok=(D>0)&(np.abs(r1/0.437-1)<0.3)&((np.abs(a/1.425-1)<0.3)|(np.abs(b/1.425-1)<0.3))&(q>=0.336)
    print(f"direcoes V0 aprovadas (q>=0,336): {int(ok.sum())}/4000 (pos {int(ok[:2000].sum())}, neg {int(ok[2000:].sum())})")
    if ok.any():
        i=np.argmax(np.where(ok,q,-1)); d=DIRS[i]; Mh={h:d@Av[h] for h in Av}; u=176.1/D[i]
        print(f"  melhor ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}) q={q[i]:.2f}; (S0-L)/D={r1[i]:.2f} (0,437); Sp1 {a[i]:.2f} Sp2 {b[i]:.2f} (1,425)")
        print(f"  massas previstas (MeV, relativas a n): Lambda {176.1:.0f}, Sigma0 {(Mh['S0']-Mh['n'])*u:.0f} (obs 253), Sigma- {(Mh['Sm']-Mh['n'])*u:.0f} (obs 258), Sigma+(1,3,5) {(Mh['Sp1']-Mh['n'])*u:.0f}, Sigma+(2,3,4) {(Mh['Sp2']-Mh['n'])*u:.0f} (obs 250), p {(Mh['p']-Mh['n'])*u:.0f} (obs -1)")
    rs=np.random.default_rng(45); okk=0;okr=0;n=40
    for _ in range(n):
        y=pert(z,rs.uniform(-1,1,4)); o=full(y)
        if o is not None and o[0]<1: okr+=1
        okk+= (f(y)<1)
    print(f"estabilidade +-10% (40): raio viavel {okr}/{n}; V0 completo (L_c<1) {okk}/{n}")
