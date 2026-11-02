# ADENDO POS-HOC (nao pre-registrado antes de ver os viaveis de V2): massas previstas, ordem do Sigma0, vacuidade.
import io, contextlib, runpy, numpy as np
with contextlib.redirect_stdout(io.StringIO()): H=runpy.run_path('P14HYB_calc.py',run_name='lib')
full=H['full']; A=H['A']; DIRS,MX=A['DIRS'],A['MX']; H_=A['H_']
for z in ([0.258,0.224,0.047,0.229],[0.080,0.439,0.0,0.083],[0.025,0.574,0.0,0.0]):
    z=np.array(z); Lr,Av,R=full(z)
    M={h:DIRS@Av[h] for h in H_}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        a=(M["Sp2"]-M["p"])/D; r1=(M["S0"]-M["L"])/D; q=D/MX
        ok=(D>0)&(np.abs(a/1.425-1)<0.3)&(q>=0.336)
    print(f"ponto {z}: L_r={Lr:.3f} R_L={R[2]:+.3f} R_Sp={R[3]:+.3f}; direcoes V2 aprovadas {int(ok.sum())}/4000 (neg {int(ok[2000:].sum())}, pos {int(ok[:2000].sum())})")
    if ok.any():
        i=np.argmax(np.where(ok,q,-1)); d=DIRS[i]
        print(f"   melhor ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}) q={q[i]:.2f}; Sigma0-Lambda prevista {r1[i]*176.1:+.0f} MeV (obs +77); (Sp-p)/D={a[i]:.2f}; fracao aprovadas com Sigma0>Lambda: {np.mean(r1[ok]>0):.2f}")
        # n-p : constante c fixa; deltas de energia estatica-dinamica
        print(f"   <A>: p {Av['p'][0]:.2f}, n {Av['n'][0]:.2f}, S0 {Av['S0'][0]:.2f}, L {Av['L'][0]:.2f}, Sp {Av['Sp2'][0]:.2f}; <B>: n {Av['n'][1]:.3f} L {Av['L'][1]:.3f}; <C>: n {Av['n'][2]:.3f} L {Av['L'][2]:.3f}")
    rng=np.random.default_rng(44); hit=0;N=300
    for _ in range(N):
        t2=rng.uniform(0.5,3.0)
        with np.errstate(all='ignore'): okr=(D>0)&(np.abs(a/t2-1)<0.3)&(q>=0.336)
        hit+=bool(okr.any())
    print(f"   vacuidade: {100*hit/N:.0f}% de alvos (Sp-p)/D aleatorios tem direcao aprovada")
