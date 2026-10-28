"""Nota (nao e teste): que Delta t cada setor precisa para o modelo geometria x escala fazer sentido."""
import math
hbar_c=197.327e-3  # GeV fm
tP=5.391e-44; lP=1.616e-35; c=2.998e8
LAM=1.2609e12; be=325/1152; S=212.677; g2=2*math.pi**2/3; K=1.16335
geo_p=0.7071/1.212**2               # fator geometrico do proton no hexagono (r_p^2 = geo_p a^2)
def dt_from_a(a_m): return a_m/c
def ell_from_a(a_m): E=hbar_c/(a_m*1e15); return be*math.log(LAM/E), E
rows=[]
# barion: a = 1,212 fm (fixado pelo proton)
a_B=1.212e-15; lB,EB=ell_from_a(a_B); rows.append(("barions (hexagono, fixado pelo proton)",a_B,lB,EB))
# eletron: compositividade > 1e5 GeV (referencia de WP-SIZE)
a_e=hbar_c/1e5*1e-15; le,Ee=ell_from_a(a_e); rows.append(("eletron (limite compositividade 1e5 GeV)",a_e,le,Ee))
# neutrino: |r2| < 1e-36 m^2 ; pior fator geometrico SZ-T: nu_e 5,4e-31 m^2 em a = 1,212 fm
fac=5.4e-31/a_B**2; a_n=math.sqrt(1e-36/fac); ln_,En=ell_from_a(a_n); rows.append(("neutrinos (limite raio de carga 1e-32 cm^2)",a_n,ln_,En))
a_L=hbar_c/LAM*1e-15; rows.append(("escala UV 1/Lambda_meta (identidade, l = 0)",a_L,0.0,LAM))
for nm,a,l,E in rows:
    print(f"{nm:46s} a {'<=' if 'limite' in nm else '='} {a:.2e} m | Delta t {'<=' if 'limite' in nm else '='} {dt_from_a(a):.2e} s = {dt_from_a(a)/tP:.2e} t_P | hbar/Delta t {'>=' if 'limite' in nm else '='} {E:.3g} GeV | degrau l {'<=' if 'limite' in nm else '='} {l:.2f}")
for l in (8,10):
    R=lP*math.exp(S/(2*g2)+l/be)/math.sqrt(K); print(f"cascata R({l}) = {R:.2e} m -> r_p previsto se a = R({l}): {math.sqrt(geo_p)*R*1e15:.3f} fm (medido 0,841)")
