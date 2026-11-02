# ADENDO POS-HOC (nao pre-registrado antes de ver os pontos viaveis conjuntos de WP-RIG): robustez, razoes exatas e naturalidade de c.
import numpy as np, runpy, io, contextlib
with contextlib.redirect_stdout(io.StringIO()): X=runpy.run_path('P14RIG_calc.py',run_name='lib')
gse,H_,DIRS,mass_ok,evalrig=X['gse'],X['H_'],X['DIRS'],X['mass_ok'],X['evalrig']
z0=np.array([0.001,0.071,0.22,0.10])
r={h:gse(h,*z0) for h in H_}; Av={h:np.array(r[h][2:5]) for h in H_}
M={h:DIRS@Av[h] for h in Av}; m=mass_ok(M)
print(f"direcoes amostradas que passam em massa: {int(m.sum())}/{len(m)}; no octante +: {int(m[:2000].sum())}; no octante -: {int(m[2000:].sum())}")
Dd=DIRS[m]; print("faixa das direcoes aprovadas (min..max por componente):",np.round(Dd.min(0),2),np.round(Dd.max(0),2))
# melhor direcao: razoes mais proximas
best=None
for d in Dd:
    Mh={h:d@Av[h] for h in H_}; D=Mh["L"]-Mh["n"]; r1=(Mh["S0"]-Mh["L"])/D; a=(Mh["Sp1"]-Mh["p"])/D; b=(Mh["Sp2"]-Mh["p"])/D
    sp=min((a,b),key=lambda v:abs(v/1.425-1)); sc=max(abs(r1/0.437-1),abs(sp/1.425-1))
    if best is None or sc<best[0]: best=(sc,d,D,r1,sp,Mh)
sc,d,D,r1,sp,Mh=best; u=176.1/D; c=1.29/u-(Mh["n"]-Mh["p"])
print(f"melhor direcao ({d[0]:+.2f},{d[1]:+.2f},{d[2]:+.2f}): D={D:.3f}, (S0-L)/D={r1:.3f} (alvo 0,437), (Sp-p)/D={sp:.3f} (alvo 1,425), desvio max {100*sc:.0f}%; |c|/D={abs(c)/D:.2f}")
print(f"   massas previstas (u={u:.1f} MeV/un): L-n={D*u:.0f}, S0-L={(Mh['S0']-Mh['L'])*u:.0f} (77), Sp-p={(Mh['Sp1' if abs((Mh['Sp1']-Mh['p'])/D/1.425-1)<abs((Mh['Sp2']-Mh['p'])/D/1.425-1) else 'Sp2']-Mh['p'])*u:.0f} (251)")
rng=np.random.default_rng(28); ok=0;okr=0;N=60
for _ in range(N):
    y=z0*(1+rng.uniform(-0.1,0.1,4)); Lc,inf=evalrig(y)
    if inf is not None and inf[2]<1: okr+=1
    if Lc<1: ok+=1
print(f"estabilidade +-10% nos 4 parametros dinamicos: raio viavel {okr}/{N}; raio viavel E massa-natural {ok}/{N}")
