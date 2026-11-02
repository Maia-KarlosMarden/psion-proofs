import math, itertools, numpy as np
Q9=np.array([1,-5,-5,10,4,4])/9
def R2(t,Q=Q9):
    th=2*math.pi*np.array(t)/6; q=Q[list(t)]; X=np.stack([np.cos(th),np.sin(th)],1); c=X.mean(0)
    v=(q*((X-c)**2).sum(1)).sum(); Qt=q.sum(); return v if abs(Qt)<1e-12 else v/Qt
rp=R2((0,4,5)); Rs=R2((0,1,2))/rp
rp2=0.8414**2; sig_fm2=math.hypot(0.12,0.09); Robs=0.61/rp2; sR=sig_fm2/rp2
print(f"SELEX (arXiv hep-ex/0106053 abstract): <r2> = 0.61 +- 0.12 (stat) +- 0.09 (syst) fm^2 -> sigma_total = {sig_fm2:.3f} fm^2 (quadrature), not 0.21 (linear sum used in the secondary source)")
print(f"K1 corrigido: R_obs = {Robs:.3f} +- {sR:.3f}; R_pred = {Rs:+.3f}; desvio = {abs(Rs-Robs)/sR:.2f} sigma; limiares 1sigma={sR:.3f}, 2sigma={2*sR:.3f}")
ok=tot=0
for perm in itertools.permutations(range(6)):
    Q=Q9[list(perm)]
    try: rr=R2((0,1,2),Q)/R2((0,4,5),Q)
    except ZeroDivisionError: continue
    if np.isfinite(rr): tot+=1; ok+=abs(rr-Robs)<=sR
print(f"acaso: fracao das {tot} permutacoes de Q9 com R_Sigma- dentro de 1 sigma corrigido: {ok/tot:.3f}")
# K2 sensibilidade
exec("M="+open('P14KILL_calc.py').read().split("M=")[1].split("\nlg=")[0])
lg=math.log(34.6); base=np.array([10-math.log(m/0.5089)/lg for m in M.values()]); fr=lambda l: np.abs(l-np.round(l))
rng=np.random.default_rng(91); out=[]
for amp in (0.01,0.05,0.20):
    ps=[]
    for _ in range(300):
        l=base-np.log(1+amp*rng.uniform(-1,1,len(base)))/lg
        Nin=(fr(l)<=0.15).sum(); off=rng.random(2000); Nn=np.array([(fr(l+o)<=0.15).sum() for o in off]); ps.append((Nn>=Nin).mean())
    print(f"K2 sensibilidade: erro aleatorio de ate +-{amp*100:.0f}% em cada massa: p mediano {np.median(ps):.2f}, minimo {np.min(ps):.2f}, maximo {np.max(ps):.2f}")
