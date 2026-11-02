import math, numpy as np
from scipy.special import spherical_jn
from scipy.optimize import brentq
from scipy.integrate import quad
hc=197.327; rp=0.841; a=1.212; Mp=938.272; LAM=1.2609e12; be=325/1152
x=brentq(lambda t: spherical_jn(0,t)-spherical_jn(1,t),1.5,3.0)
rho=lambda u: spherical_jn(0,u)**2+spherical_jn(1,u)**2   # u = k r, k = x/R
num=quad(lambda u: u**4*rho(u),0,x)[0]; den=quad(lambda u: u**2*rho(u),0,x)[0]
r2R2=num/den/x**2            # <r^2>/R^2
print(f"x = {x:.4f}; <r^2>/R^2 = {r2R2:.4f} (rms = {math.sqrt(r2R2):.4f} R)")
R3=rp/math.sqrt(r2R2)
rows=[("A0 E=hbar c/a",hc/a),("A1 E=hbar c/(2a)",hc/(2*a)),("A2 bag, R=a",x*hc/a),(f"A3 bag, R=r_p/rms ({R3:.3f} fm)",x*hc/R3)]
print("O1: M_p = 3E vs 938,272 MeV")
for nm,E in rows:
    M=3*E; print(f"   {nm:32s}: E = {E:7.1f} MeV, M_p = {M:7.1f} MeV ({M/Mp-1:+.0%}) -> {'passa' if abs(M/Mp-1)<=0.10 else 'falha'}")
print("O2: orbital comum para os tres psions: R_n = 0, R_Sigma- = 1 -> vs medidos -0,164 e 0,86+-0,21: R_n falha (0 vs -0,164 e -100% de desvio), R_S passa")
print("O3 (informativo): E por psion vs degraus:")
for l in (8,8.5):
    m=LAM*math.exp(-l/be); print(f"   l={l}: {m*1e3:.0f} MeV; razao aos E de A0..A3: "+", ".join(f"{m*1e3/E:.2f}" for _,E in rows))
# look-elsewhere: janela 10% em fator log-uniforme [0.3,3]
w=2*math.log(1.1)/math.log(10); print(f"look-elsewhere: janela = {w:.3f} do intervalo; 4 correcoes: {1-(1-w)**4:.2f}")
