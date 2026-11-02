# ADENDO POS-HOC (nao pre-registrado antes de ver os otimos robustos de HYB3): massas, direcao natural, Sigma- a parte.
import io, contextlib, runpy, numpy as np
with contextlib.redirect_stdout(io.StringIO()): H=runpy.run_path('P14HYB3_calc.py',run_name='lib')
full3,A,static,gse=H['full3'],H['A'],H['static'],H['gse']; DIRS,MX=A['DIRS'],A['MX']
H2=runpy.run_path('P14HYB2_calc.py',run_name='lib'); DY=H2['DY']
for z in ([0.258,0.224,0.047,0.229],[0.080,0.439,0.0,0.083],[0.035,0.456,0.095,0.0]):
    z=np.array(z); R,Av=full3(z)
    M={h:DIRS@Av[h] for h in Av}; D=M["L"]-M["n"]
    with np.errstate(all='ignore'):
        r1=(M["S0"]-M["L"])/D; a=(M["Sp1"]-M["p"])/D; b=(M["Sp2"]-M["p"])/D; q=D/MX
        ok=(D>0)&(np.abs(r1/0.437-1)<0.3)&((np.abs(a/1.425-1)<0.3)|(np.abs(b/1.425-1)<0.3))&(q>=0.336)
    print(f"ponto {z}: R_n={R[0]:+.3f} R_S0={R[1]:+.3f} R_L={R[2]:+.3f} R_Sp1={R[3]:+.3f} R_Sp2={R[4]:+.3f}; direcoes V0 aprovadas {int(ok.sum())}/4000 (pos {int(ok[:2000].sum())}, neg {int(ok[2000:].sum())})")
    if ok.any():
        i=np.argmax(np.where(ok,q,-1)); d=DIRS[i]; Mh={h:d@Av[h] for h in Av}; u=176.1/D[i]
        print(f"  melhor ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}) q={q[i]:.2f}; (S0-L)/D={r1[i]:.2f}; Sp1 {a[i]:.2f} Sp2 {b[i]:.2f}; faixa (S0-L)/D nas aprovadas [{np.nanmin(r1[ok]):.2f}, {np.nanmax(r1[ok]):.2f}]")
        print(f"  massas rel. a n (MeV): Lambda 176, Sigma0 {(Mh['S0']-Mh['n'])*u:.0f} (253), Sigma+(2,3,4) {(Mh['Sp2']-Mh['n'])*u:.0f} (250), Sigma+(1,3,5) {(Mh['Sp1']-Mh['n'])*u:.0f}, p {(Mh['p']-Mh['n'])*u:+.0f}")
        # Sigma- a parte
        t,s,k,l=z; sm=gse("Sm",t,s,k,l); rp=static("p")[0]; rs,As=static("Sm"); Asm=np.array(sm[2:5])
        print(f"  Sigma- (fora do laco): dinamico R={sm[0]/rp:+.2f}, massa rel. n {(d@Asm-Mh['n'])*u:.0f} MeV; estatico R={rs/rp:+.2f}, massa rel. n {(d@As-Mh['n'])*u:.0f} MeV (obs R 0,86+-0,21; massa 258)")
