# ADENDO POS-HOC (nao pre-registrado antes de ver os viaveis de V2/V3/V4/V5 no DYN12): detalhes e estabilidade do ponto de V2.
import io, contextlib, runpy, numpy as np
with contextlib.redirect_stdout(io.StringIO()): A=runpy.run_path('P14ABL_calc.py',run_name='lib')
full_dy,DIRS,MX,H_,qvar,SPEC=A['full_dy'],A['DIRS'],A['MX'],A['H_'],A['qvar'],A['SPEC']
for tag,z in (("V2",[0.478,0.156,0.0,0.48]),("V4",[0.676,0.039,0.0,0.595])):
    z=np.array(z); Lr,Av=full_dy(z); g=A['DY']['gse']
    r={h:g(h,*z) for h in H_}; rp=r["p"][0]; R=[r[h][0]/rp for h in ("n","S0","L","Sp1","Sp2")]
    print(f"{tag} ponto {z}: L_r={Lr:.3f}; R_n={R[0]:+.3f} R_S0={R[1]:+.3f} R_L={R[2]:+.3f} R_Sp1={R[3]:+.3f} R_Sp2={R[4]:+.3f}")
    M={h:DIRS@Av[h] for h in H_}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        r2=np.minimum(np.abs((M["Sp1"]-M["p"])/D/1.425-1),np.abs((M["Sp2"]-M["p"])/D/1.425-1)); r1=(M["S0"]-M["L"])/D
        ok=(D>0)&(r2<0.3)&(D/MX>=0.336)
    print(f"  direcoes que passam em V2 com q>=0,336: {int(ok.sum())}/4000; (S0-L)/D nessas: mediana {np.median(r1[ok]):.2f} min {r1[ok].min():.2f} max {r1[ok].max():.2f} (alvo 0,437)" if ok.any() else "  nenhuma direcao")
    if ok.any():
        i=np.argmax(np.where(ok,D/MX,-1)); d=DIRS[i]; print(f"  melhor direcao ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}) q={D[i]/MX[i]:.3f}; (S0-L)/D={r1[i]:.2f}")
    rng=np.random.default_rng(37); okc=0;okr=0;N=60
    for _ in range(N):
        y=np.clip(z*(1+rng.uniform(-0.1,0.1,4)),0,None); out=full_dy(y)
        if out is None: continue
        L2,Av2=out; okr+= (L2<1); q2=qvar(Av2,SPEC["V2"]); okc+= (L2<1 and q2 is not None and q2>=0.336)
    print(f"  estabilidade +-10%: raio viavel {okr}/{N}; raio e V2 (q>=0,336) {okc}/{N}")
