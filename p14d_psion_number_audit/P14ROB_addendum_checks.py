# ADENDO POS-HOC (nao pre-registrado antes de ver o viavel robusto): massas previstas, direcoes, vacuidade do criterio V2 no ponto robusto.
import io, contextlib, runpy, numpy as np
with contextlib.redirect_stdout(io.StringIO()): A=runpy.run_path('P14ABL_calc.py',run_name='lib')
full_dy,DIRS,MX,H_=A['full_dy'],A['DIRS'],A['MX'],A['H_']; g=A['DY']['gse']
z=np.array([0.001,0.065,0.287,0.095]); Lr,Av=full_dy(z)
r={h:g(h,*z) for h in H_}; rp=r["p"][0]; R=[r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]
print(f"ponto {z}: L_r={Lr:.3f} R_n={R[0]:+.3f} R_S0={R[1]:+.3f} R_L={R[2]:+.3f} R_Sp1={R[3]:+.3f} R_Sp2={R[4]:+.3f}")
M={h:DIRS@Av[h] for h in H_}; D=M["L"]-M["n"]
with np.errstate(all='ignore'):
    a=(M["Sp1"]-M["p"])/D; b=(M["Sp2"]-M["p"])/D; r1=(M["S0"]-M["L"])/D; q=D/MX
    ok=(D>0)&((np.abs(a/1.425-1)<0.3)|(np.abs(b/1.425-1)<0.3))&(q>=0.336)
print(f"direcoes que passam em V2 com q>=0,336: {int(ok.sum())}/4000 (positivas {int(ok[:2000].sum())}, negativas {int(ok[2000:].sum())})")
if ok.any():
    i=np.argmax(np.where(ok,q,-1)); d=DIRS[i]; u=176.1/D[i]
    print(f"  melhor direcao ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}) q={q[i]:.3f}; (S0-L)/D={r1[i]:+.2f} => Sigma0-Lambda prevista {r1[i]*176.1:+.0f} MeV (obs +77); (Sp-p)/D: Sp1 {a[i]:.2f}, Sp2 {b[i]:.2f} (alvo 1,425)")
    print(f"  (S0-L)/D entre as direcoes aprovadas: min {np.nanmin(r1[ok]):+.2f}, mediana {np.nanmedian(r1[ok]):+.2f}, max {np.nanmax(r1[ok]):+.2f}")
    print(f"  fracao das aprovadas com Sigma0>Lambda (ordem correta): {np.mean(r1[ok]>0):.2f}")
# vacuidade: para alvos r2 aleatorios, existe direcao natural com q>=0,336 neste ponto?
rng=np.random.default_rng(42); hit=0;N=300
for _ in range(N):
    t2=rng.uniform(0.5,3.0)
    with np.errstate(all='ignore'):
        okr=(D>0)&((np.abs(a/t2-1)<0.3)|(np.abs(b/t2-1)<0.3))&(q>=0.336)
    hit+=bool(okr.any())
print(f"vacuidade: alvos (Sp-p)/D aleatorios U[0,5;3,0]: {hit}/{N} ({100*hit/N:.0f}%) tem alguma direcao natural com q>=0,336 neste ponto (alvo real: {'sim' if ok.any() else 'nao'})")
# faixa de r2 alcancavel
with np.errstate(all='ignore'):
    m=(D>0)&(q>=0.336); vals=np.concatenate([a[m],b[m]]); vals=vals[np.isfinite(vals)]
print(f"faixa de (Sp-p)/D alcancavel com D>0 e q>=0,336: [{vals.min():.2f}, {vals.max():.2f}]")
