import math
from fractions import Fraction as F
hc=0.197327; LAM=1.2609e12; be=325/1152; geo=0.7071/1.212**2
Q={k:F(q,9) for k,q in enumerate((1,-5,-5,10,4,4))}
print(f"A1 (geo_p = {geo:.4f}):")
best=None
for l2 in range(0,41):
    l=l2/2; m=LAM*math.exp(-l/be); rp=math.sqrt(geo)*hc/m
    if 0.05<rp<20: print(f"   l = {l:4.1f}: m(l) = {m:.3e} GeV, r_p = {rp:.3f} fm")
    if best is None or abs(math.log(rp/0.841))<abs(math.log(best[1]/0.841)): best=(l,rp)
print(f"   melhor: l = {best[0]} -> r_p = {best[1]:.3f} fm (desvio {best[1]/0.841-1:+.0%}) -> {'passa' if abs(best[1]/0.841-1)<0.1 else 'falha'}")
p,n=(0,4,5),(0,1,5)
q2=lambda t: sum(Q[k]**2 for k in t)
SG=q2(p)-q2(n); SGq=F(4,9)*2+F(1,9)-(F(4,9)+F(1,9)*2)
print(f"A2 Gottfried valencia: PH S_G = {q2(p)} - {q2(n)} = {SG} = {float(SG):+.3f}; quarks = {SGq} = {float(SGq):.3f}; medido 0,235 +- 0,026")
ok=(SG>0) and 0.5<float(SG)/0.235<2
print(f"   -> {'passa' if ok else 'falha'}; mar assimetrico exigido para chegar ao medido: {0.235-float(SG):+.3f} (quarks: {0.235-float(SGq):+.3f})")
print(f"A3 controle: n de partons = 3 e 3; carga p = {sum(Q[k] for k in p)}, n = {sum(Q[k] for k in n)}")
# B
QB=LAM**2; print(f"B1 escala de queda prevista Q^2 ~ Lambda_meta^2 = {QB:.2e} GeV^2 vs 0,71 GeV^2: fator {QB/0.71:.1e} -> falha")
print(f"B2 resolucao dos psions em Q ~ {LAM:.1e} GeV vs partons observados em 1-100 GeV -> falha")
